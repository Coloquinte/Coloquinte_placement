
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
        for(index_t l=c.neighbours_begin; l<c.neighbours_begin + c.height; ++l){
            // not verified now since we don't modify the position for the obstacles
            // : assert(c.position.x_ >= min_x_ and c.position.x_ + c.width <= max_x_);

            assert(c.row + c.height <= row_cnt());

            if(neighbours_[l].first != null_ind){
                auto const oc = cells_[neighbours_[l].first];
                // Correct neighbour position
                assert(c.position.x_ >= oc.position.x_ + oc.width);
            }
            else{
                // Beginning of a row
                //assert(row_first_cells_[c.row + l - c.neighbours_begin] == i);
                assert(row_first_cells_[c.row + l - c.neighbours_begin] == i);
            }
            if(neighbours_[l].second != null_ind){
                auto const oc = cells_[neighbours_[l].second];
                // Correct neighbour position
                assert(c.position.x_ + c.width <= oc.position.x_);
            }
            else{
                // End of a row
                assert(row_last_cells_[c.row + l - c.neighbours_begin] == i);
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

} // namespace dp
} // namespace coloquinte


