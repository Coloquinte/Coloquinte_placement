
#include "Coloquinte/detailed.hxx"

#include <lemon/smart_graph.h>
#include <lemon/network_simplex.h>

#include <cassert>

namespace coloquinte{
namespace dp{

detailed_placement::detailed_placement(
        std::vector<internal_cell> const cells,
        std::vector<std::vector<index_t> > const rows,
        int_t min_x, int_t max_x,
        int_t y_origin,
        index_t nbr_rows, int_t row_height
    )
    :
        min_x_(min_x), max_x_(max_x),
        y_origin_(y_origin),
        cells_(cells)
    {

    assert(row_height > 0);
    assert(min_x < max_x);
    assert(rows.size() == nbr_rows);

    index_t nbr_lims = 0; 
    for(internal_cell & c : cells_){
        c.neighbours_begin = nbr_lims;
        nbr_lims += c.height;
    }

    neighbours_ .resize(nbr_lims, std::pair<index_t, index_t>(null_ind, null_ind) );

    row_first_cells_ .resize(nbr_rows, null_ind);
    row_last_cells_  .resize(nbr_rows, null_ind);

    std::vector<bool> explored(nbr_lims, false);
    // Now we extract the dependencies
    for(index_t r=0; r<rows.size(); ++r){

        if(not rows[r].empty()){
            row_first_cells_[r] = rows[r].front();
            row_last_cells_[r]  = rows[r].back();
        }

        for(index_t c : rows[r]){
            // Has this row of the cell already been visited?
            assert(not explored[neighbour_index(c, r)]);
            explored[neighbour_index(c, r)] = true;
        }

        for(index_t i=0; i+1<rows[r].size(); ++i){
            index_t c1 = rows[r][i], c2 = rows[r][i+1];

            // Save in the internal format
            neighbours_[neighbour_index(c1, r)].second = c2;
            neighbours_[neighbour_index(c2, r)].first  = c1;

            // The positions are correct
            assert(cells_[c1].position.x_ + cells_[c1].width <= cells_[c2].position.x_ );
        }
    }

    // Every level of every cell must have been visited
    for(bool o : explored)
        assert(o);

    // Verify that we haven't made any obvious mistake
    selfcheck();
}

void detailed_placement::selfcheck() const{
    assert(row_first_cells_.size() == row_last_cells_.size());

    for(index_t i=0; i<cell_cnt(); ++i){
        internal_cell const c = cells_[i];
        for(index_t l=0; l<c.height; ++l){
            // not verified now since we don't modify the position for the obstacles
            // : assert(c.position.x_ >= min_x_ and c.position.x_ + c.width <= max_x_);

            index_t n_ind = l + c.neighbours_begin;
            assert(c.row + c.height <= row_cnt());

            if(neighbours_[n_ind].first != null_ind){
                index_t oi = neighbours_[n_ind].first;
                auto const oc = cells_[oi];
                // Correct neighbour position
                assert(c.position.x_ >= oc.position.x_ + oc.width);
                assert(neighbours_[neighbour_index(oi, c.row+l)].second == i);
            }
            else{
                // Beginning of a row
                //assert(row_first_cells_[c.row + l - c.neighbours_begin] == i);
                assert(row_first_cells_[c.row + n_ind - c.neighbours_begin] == i);
            }
            if(neighbours_[n_ind].second != null_ind){
                index_t oi = neighbours_[n_ind].second;
                auto const oc = cells_[oi];
                // Correct neighbour position
                assert(c.position.x_ + c.width <= oc.position.x_);
                assert(neighbours_[neighbour_index(oi, c.row+l)].first == i);
            }
            else{
                // End of a row
                assert(row_last_cells_[c.row + n_ind - c.neighbours_begin] == i);
            }
        }
    }
}

void optimize_positions(netlist const & circuit, detailed_placement & pl){
    // Solves a minimum cost flow problem to optimize the placement at fixed topology
    // Concretely, it means aligning the pins to minimize the wirelength
    // It uses the Lemon network simplex solver from the Coin-OR initiative, which should scale well up to hundred of thousands of cells

    using namespace lemon;
    DIGRAPH_TYPEDEFS(SmartDigraph);
    // Create a graph with the cells and bounds of the nets as node
    SmartDigraph g;

    std::vector<Node> cell_nodes(circuit.cell_cnt());
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        if((circuit.get_cell(i).attributes & XMovable) != 0)
            cell_nodes[i] = g.addNode();
    }
    std::vector<Node> Lnet_nodes(circuit.net_cnt()), Unet_nodes(circuit.net_cnt());
    for(index_t i=0; i<circuit.net_cnt(); ++i){
        Lnet_nodes[i] = g.addNode();
        Unet_nodes[i] = g.addNode();
    }

    // Two nodes for position constraints
    Node fixed = g.addNode();

    typedef std::pair<SmartDigraph::Arc, int_t> arc_pair;
    typedef std::pair<SmartDigraph::Node, int_t> node_pair;
    // The arcs corresponding to constraints of the original problem
    std::vector<arc_pair> constraint_arcs;

    // Now we add every positional constraint, which becomes an arc in the min-cost flow problem
    for(index_t i=0; i<circuit.cell_cnt(); ++i){ // The cells
        auto c = pl.cells_[i];
        for(index_t l = c.neighbours_begin; l < c.neighbours_begin + c.height; ++l){
            index_t oi = pl.neighbours_[l].second;
            if(oi == null_ind) continue;
            auto oc = pl.cells_[oi];
            assert(c.position.x_ + c.width <= oc.position.x_);
            
            if((circuit.get_cell(i).attributes & XMovable) != 0 and (circuit.get_cell(oi).attributes & XMovable) != 0){
                // Two movable cells: OK
                auto A = g.addArc(cell_nodes[oi], cell_nodes[i]);
                constraint_arcs.push_back(arc_pair(A, -c.width));
            }
            else if((circuit.get_cell( i).attributes & XMovable) != 0){
                // The cell c is movable and constrained on the right
                auto A = g.addArc(fixed, cell_nodes[i]);
                constraint_arcs.push_back(arc_pair(A, oc.position.x_ - c.width));
            }
            else if((circuit.get_cell(oi).attributes & XMovable) != 0){
                // The cell oc is movable and constrained on the left
                auto A = g.addArc(cell_nodes[oi], fixed);
                constraint_arcs.push_back(arc_pair(A, -c.position.x_ - c.width));
            }
        }
    }

    
    for(index_t r=0; r<pl.row_cnt(); ++r){ // And the boundaries of each row
        index_t lc = pl.row_first_cells_[r];
        if(lc != null_ind and (circuit.get_cell(lc).attributes & XMovable) != 0){
            auto Al = g.addArc(cell_nodes[lc], fixed);
            constraint_arcs.push_back(arc_pair(Al, -pl.min_x_));
        }
    }
    for(index_t r=0; r<pl.row_cnt(); ++r){ // And the boundaries of each row
        index_t rc = pl.row_last_cells_[r];
        if(rc != null_ind and (circuit.get_cell(rc).attributes & XMovable) != 0){
            auto Ar = g.addArc(fixed, cell_nodes[rc]);
            constraint_arcs.push_back(arc_pair(Ar, pl.max_x_ - pl.cells_[rc].width));
        }
    }
    

    // And every pin of every net: arcs too
    for(index_t n=0; n<circuit.net_cnt(); ++n){
        assert(circuit.get_net(n).pin_cnt > 0);
        for(auto p : circuit.get_net(n)){
            index_t c = p.cell_ind;
            int_t pin_offs = 0.5 * pl.cells_[c].width + (pl.cells_[c].x_orientation ? p.offset.x_ : - p.offset.x_); // Offset to the beginning of the cell
            if((circuit.get_cell(c).attributes & XMovable) != 0){
                Arc Al = g.addArc(cell_nodes[c], Lnet_nodes[n]);
                constraint_arcs.push_back(arc_pair(Al, pin_offs));
                Arc Ar = g.addArc(Unet_nodes[n], cell_nodes[c]);
                constraint_arcs.push_back(arc_pair(Ar, -pin_offs));
            }
            else{ // Fixed offset
                auto Al = g.addArc(fixed, Lnet_nodes[n]);
                constraint_arcs.push_back(arc_pair(Al, pl.cells_[c].position.x_ + pin_offs));
                auto Ar = g.addArc(Unet_nodes[n], fixed);
                constraint_arcs.push_back(arc_pair(Ar, - pl.cells_[c].position.x_ - pin_offs));
            }
        }
    }

    // Then the only capacitated arcs: the ones for the nets
    std::vector<node_pair> net_supplies;
    for(index_t n=0; n<circuit.net_cnt(); ++n){
        // TODO: use net weights
        net_supplies.push_back(node_pair(Unet_nodes[n], 1));
        net_supplies.push_back(node_pair(Lnet_nodes[n], -1));
    }

    // Create the maps to have cost and capacity for the arcs
    IntArcMap cost(g, 0);
    IntArcMap capacity(g, circuit.net_cnt());
    IntNodeMap supply(g, 0);

    for(arc_pair A : constraint_arcs){
        cost[A.first] = A.second;
    }

    for(node_pair N : net_supplies){
        supply[N.first] = N.second;
    }

    // Then we (hope the solver can) solve it
    NetworkSimplex<SmartDigraph> ns(g);
    ns.supplyMap(supply).costMap(cost);
    auto res = ns.run();
    if(res != ns.OPTIMAL){
        abort();
    }
    
    // And we get the new positions as the dual values of the current solution (compared to the fixed pin) 
    for(index_t c=0; c<circuit.cell_cnt(); ++c){ // The cells
        if((circuit.get_cell(c).attributes & XMovable) != 0){
            pl.cells_[c].position.x_ = ns.potential(cell_nodes[c]) - ns.potential(fixed);
        }
    }
    pl.selfcheck();
}

namespace{

inline float_t get_nets_cost(netlist const & circuit, detailed_placement const & pl, std::vector<index_t> const & involved_nets){
    float_t cost = 0.0;

    for(index_t n : involved_nets){
        if(circuit.get_net(n).pin_cnt == 0) continue;
        float const INF = std::numeric_limits<float_t>::infinity();
        auto mins = point<float_t>( INF,  INF);
        auto maxs = point<float_t>(-INF, -INF);
        for(netlist::pin_t p : circuit.get_net(n)){
            auto offs = point<float_t>(
                pl.cells_[p.cell_ind].x_orientation ? p.offset.x_ : -p.offset.x_,
                pl.cells_[p.cell_ind].y_orientation ? p.offset.y_ : -p.offset.y_
            );
            point<float_t> pos = offs + static_cast<point<float_t> >(pl.cells_[p.cell_ind].position) + static_cast<float_t>(0.5) * static_cast<point<float_t> >(circuit.get_cell(p.cell_ind).size);
            mins.x_ = std::min(mins.x_, pos.x_);
            mins.y_ = std::min(mins.y_, pos.y_);
            maxs.x_ = std::max(maxs.x_, pos.x_);
            maxs.y_ = std::max(maxs.y_, pos.y_);
        }
        cost += (maxs.x_ - mins.x_) + (maxs.y_ - mins.y_);
    }

    return cost;
}

// TODO: take fixed cells into account, even if usually one row <==> movable
// Tries to swap two cells; 
bool try_swap(netlist const & circuit, detailed_placement & pl, index_t c1, index_t c2){
    index_t row_c1 = pl.cells_[c1].row,
            row_c2 = pl.cells_[c2].row;

    assert(pl.cell_hght(c1) == 1 and pl.cell_hght(c2) == 1);
    assert(circuit.get_cell(c1).size.y_ == circuit.get_cell(c2).size.y_); // Same (standard cell) height
    assert(row_c1 != row_c2); // Same (standard cell) height
    assert( (circuit.get_cell(c1).attributes & XMovable) != 0 and (circuit.get_cell(c1).attributes & YMovable) != 0);
    assert( (circuit.get_cell(c2).attributes & XMovable) != 0 and (circuit.get_cell(c2).attributes & YMovable) != 0);

    int_t c1_l, c1_u,
          c2_l, c2_u;

    index_t b_c1 = pl.neighbours_[pl.cells_[c1].neighbours_begin].first;
    index_t b_c2 = pl.neighbours_[pl.cells_[c2].neighbours_begin].first;
    index_t a_c1 = pl.neighbours_[pl.cells_[c1].neighbours_begin].second;
    index_t a_c2 = pl.neighbours_[pl.cells_[c2].neighbours_begin].second;

    // Limit positions for c1 and c2
    c1_l = b_c1 != null_ind ? pl.cells_[b_c1].position.x_ + pl.cells_[b_c1].width : pl.min_x_;
    c2_l = b_c2 != null_ind ? pl.cells_[b_c2].position.x_ + pl.cells_[b_c2].width : pl.min_x_;
    c1_u = a_c1 != null_ind ? pl.cells_[a_c1].position.x_ : pl.max_x_;
    c2_u = a_c2 != null_ind ? pl.cells_[a_c2].position.x_ : pl.max_x_;

    int_t swp_min_c1, swp_max_c1,
          swp_min_c2, swp_max_c2;

    // Get the possible positions for a swap
    swp_min_c1 = c2_l;
    swp_min_c2 = c1_l;
    swp_max_c1 = c2_u - pl.cells_[c1].width;
    swp_max_c2 = c1_u - pl.cells_[c2].width;

    if(swp_max_c1 >= swp_min_c1 and swp_max_c2 >= swp_min_c2){
        // Check both orientations of the cell

        // Get all the nets involved and uniquify them (nets with more than one pin on the cells)
        std::vector<index_t> involved_nets;
        for(netlist::pin_t p : circuit.get_cell(c1)){
            involved_nets.push_back(p.net_ind);
        }
        for(netlist::pin_t p : circuit.get_cell(c2)){
            involved_nets.push_back(p.net_ind);
        }
        std::sort(involved_nets.begin(), involved_nets.end());
        involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

        // Test the cost for the old position and the cost swapping the cells
        float_t old_cost = get_nets_cost(circuit, pl, involved_nets);

        // Save the old values
        int_t c1_x = pl.cells_[c1].position.x_;
        int_t c2_x = pl.cells_[c2].position.x_;
        int_t c1_y = pl.cells_[c1].position.y_;
        int_t c2_y = pl.cells_[c2].position.y_;

        // Warning: won't work if the two cells don't have the same height
        pl.cells_[c1].position.x_ = (swp_min_c1 + swp_max_c1) / 2;
        pl.cells_[c2].position.x_ = (swp_min_c2 + swp_max_c2) / 2;
        pl.cells_[c1].position.y_ = c2_y;
        pl.cells_[c2].position.y_ = c1_y;

        float_t swp_cost = get_nets_cost(circuit, pl, involved_nets);

        if(swp_cost < old_cost){
            // Update the cells' neighbours
            std::swap(pl.neighbours_[pl.cells_[c1].neighbours_begin],  pl.neighbours_[pl.cells_[c2].neighbours_begin]);

            // Update the neighbours too
            if(b_c1 != null_ind) pl.neighbours_[pl.neighbour_index(b_c1, row_c1)].second = c2;
            else pl.row_first_cells_[row_c1] = c2;
            if(b_c2 != null_ind) pl.neighbours_[pl.neighbour_index(b_c2, row_c2)].second = c1;
            else pl.row_first_cells_[row_c2] = c1;
            if(a_c1 != null_ind) pl.neighbours_[pl.neighbour_index(a_c1, row_c1)].first  = c2;
            else pl.row_last_cells_[row_c1] = c2;
            if(a_c2 != null_ind) pl.neighbours_[pl.neighbour_index(a_c2, row_c2)].first  = c1;
            else pl.row_last_cells_[row_c2] = c1;

            // Update the rows
            std::swap(pl.cells_[c1].row, pl.cells_[c2].row);

            // We kept the swap
            return true;
        }
        else{
            // Reset the old values
            pl.cells_[c1].position.x_ = c1_x;
            pl.cells_[c2].position.x_ = c2_x;
            pl.cells_[c1].position.y_ = c1_y;
            pl.cells_[c2].position.y_ = c2_y;

            // We didn't swap
            return false;
        }

        // A better solution would be
        // Check the cost on y depending on the position (extremely simple: two positions for each cell)
        // Check the cost on x depending on the position: piecewise linear and relatively complex
        //      * Get all external pins
        //      * Get all nets involving only one of the cells: piecewise linear cost for each of them
        //      * For nets involving the two cells, we have an additional cost

    }
    else{ // We just cannot swap those two cells without pushing anything
        return false;
    }
}

inline index_t get_first_standard_cell(detailed_placement const & pl, index_t r, index_t c){
    while(c != null_ind and pl.cells_[c].height != 1){
        index_t next_c = pl.neighbours_[pl.neighbour_index(c, r)].second;
        assert(c != next_c);
        c = next_c;
    }
    assert(c == null_ind or pl.cells_[c].row == r);
    return c;
}

inline index_t get_first_cell_on_row(detailed_placement const & pl, index_t r){
    index_t c = pl.row_first_cells_[r];
    return get_first_standard_cell(pl, r, c);
}

inline index_t get_next_cell_on_row(detailed_placement const & pl, index_t c){
    assert(pl.cells_[c].height == 1);
    index_t next_c = pl.neighbours_[pl.cells_[c].neighbours_begin].second;
    assert(next_c != c);
    index_t r = pl.cells_[c].row;
    index_t ret = get_first_standard_cell(pl, r, next_c);
    assert(ret != c);
    return ret;
}

} // End anonymous namespace


void optimize_swaps(netlist const & circuit, detailed_placement & pl, index_t row_extent, index_t cell_extent){
    for(index_t main_row = 0; main_row < pl.row_cnt(); ++main_row){

        for(index_t other_row = main_row+1; other_row <= std::min(pl.row_cnt()-1, main_row+row_extent) ; ++other_row){

            index_t first_oc = get_first_cell_on_row(pl, other_row); // The first candidate cell to be examined
            for(index_t c = get_first_cell_on_row(pl, main_row); c != null_ind; c = get_next_cell_on_row(pl, c)){
                assert(pl.cells_[c].row == main_row);
                // Number of cells after/before the end of the cell
                index_t nb_after  = 0;
                index_t nb_before = 0;
                int_t pos_low = pl.cells_[c].position.x_ - pl.cells_[c].width, pos_hgh = pl.cells_[c].position.x_ + 2*pl.cells_[c].width;
                for(index_t oc=first_oc; oc != null_ind and nb_after <= row_extent; oc = get_next_cell_on_row(pl, oc)){
                    assert(pl.cells_[oc].row == other_row);

                    // Count the cells which should trigger stop or shouldn't be used at the next iteration
                    if(pl.cells_[oc].position.x_ >= pos_hgh) ++nb_after;
                    if(pl.cells_[oc].position.x_ + pl.cells_[oc].width <= pos_low) ++ nb_before;

                    if(try_swap(circuit, pl, c, oc)){
                        std::swap(c, oc);
                        if(c == first_oc) first_oc = oc;
                    }
                }
                while(nb_before > cell_extent){
                    nb_before--;
                    first_oc = get_next_cell_on_row(pl, first_oc);
                }
            }
        }
    }
    pl.selfcheck();
}


} // namespace dp
} // namespace coloquinte


