
#include "Coloquinte/detailed.hxx"

#include <lemon/smart_graph.h>
#include <lemon/cost_scaling.h>
#include <lemon/capacity_scaling.h>
#include <lemon/network_simplex.h>

#include <cassert>

#include <iostream>

namespace coloquinte{
namespace dp{

detailed_placement::detailed_placement(
        std::vector<int_t>    positions,
        std::vector<index_t>  row_positions,
        std::vector<int_t>    widths,
        std::vector<index_t>  heights,
        std::vector<bool> x_orientations,
        std::vector<bool> y_orientations,
        std::vector<std::vector<index_t> > const rows,
        int_t min_x, int_t max_x,
        int_t y_origin,
        index_t nbr_rows, int_t row_height
    )
    :
        min_x_(min_x), max_x_(max_x),
        y_origin_(y_origin),
        row_height_(row_height),
        positions_(positions),
        cell_rows_(row_positions),
        widths_(widths),
        x_orientations_(x_orientations),
        y_orientations_(y_orientations)
    {

    index_t sz = positions.size();

    assert(row_height > 0);
    assert(min_x < max_x);
    assert(row_positions.size() == sz
        and widths.size()  == positions.size()
        and heights.size() == positions.size()
    );
    assert(rows.size() == nbr_rows);

    cell_lims_.resize(1, 0);
    index_t nbr_lims = 0; 
    for(index_t h : heights){
        nbr_lims += h;
        cell_lims_.push_back(nbr_lims);
    }

    cells_before_ .resize(nbr_lims, null_ind);
    cells_after_  .resize(nbr_lims, null_ind);

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
            // The row_positions are correct
            assert( r >= cell_rows_[c] and r < cell_rows_[c] + cell_hght(c));

            // Has this row of the cell already been visited?
            assert(not explored[cell_lims_[c] + r - cell_rows_[c]]);
            explored[cell_lims_[c] + r - cell_rows_[c]] = true;
        }

        for(index_t i=0; i+1<rows[r].size(); ++i){
            index_t c1 = rows[r][i], c2 = rows[r][i+1];
            assert(c1 < 1000000 && c2 < 1000000);

            // Save in the internal format
            cells_after_ [cell_lims_[c1] + r - cell_rows_[c1]] = c2;
            cells_before_[cell_lims_[c2] + r - cell_rows_[c2]] = c1;

            // The positions are correct
            assert( positions[c1] + widths[c1] <= positions[c2] );
        }
    }

    // Every level of every cell must have been visited
    for(bool o : explored)
        assert(o);

    // Verify that we haven't made any obvious mistake
    selfcheck();
}

void detailed_placement::selfcheck() const{
    assert(cells_before_.size()   == cells_after_.size());
    assert( cell_rows_.size()     == cell_cnt()
        and widths_.size()        == cell_cnt()
        and positions_.size()     == cell_cnt()
        and x_orientations_.size() == cell_cnt()
        and y_orientations_.size() == cell_cnt()
    );
    assert(row_first_cells_.size() == row_last_cells_.size());

    for(index_t i=0; i<cell_cnt(); ++i){
        for(index_t l=cell_lims_[i]; l<cell_lims_[i+1]; ++l){
            // Correct x position
            assert(positions_[i] >= min_x_ and positions_[i] + widths_[i] <= max_x_);
            assert(cell_rows_[i] + cell_hght(i) <= row_cnt());

            if(cells_before_[l] != null_ind){
                // Correct neighbour position
                assert(positions_[i] >= positions_[cells_before_[l]] + widths_[cells_before_[l]]);
            }
            else{
                // Beginning of a row
                assert(row_first_cells_[cell_rows_[i] + l - cell_lims_[i]] == i);
            }
            if(cells_after_[l] != null_ind){
                // Correct neighbour position
                assert(positions_[i] + widths_[i] <= positions_[cells_after_[l]]);
            }
            else{
                // End of a row
                assert(row_last_cells_[cell_rows_[i] + l - cell_lims_[i]] == i);
            }
        }
    }
    
    for(index_t i=0; i<row_cnt(); ++i){
        if(row_first_cells_[i] != null_ind){
            assert(cells_before_[cell_lims_[row_first_cells_[i]] + i - cell_rows_[row_first_cells_[i]]] == null_ind);
        }
        if(row_last_cells_[i] != null_ind){
            assert(cells_after_[cell_lims_[row_last_cells_[i]] + i - cell_rows_[row_last_cells_[i]]] == null_ind);
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
    for(index_t c=0; c<circuit.cell_cnt(); ++c){ // The cells
        for(index_t l = pl.cell_lims_[c]; l < pl.cell_lims_[c+1]; ++l){
            index_t oc = pl.cells_after_[l];
            if(oc == null_ind) continue;

            assert(pl.positions_[c] + pl.widths_[c] <= pl.positions_[oc]);
            
            if((circuit.get_cell(c).attributes & XMovable) != 0 and (circuit.get_cell(oc).attributes & XMovable) != 0){
                // Two movable cells: OK
                auto A = g.addArc(cell_nodes[oc], cell_nodes[c]);
                constraint_arcs.push_back(arc_pair(A, -pl.widths_[c]));
            }
            else if((circuit.get_cell( c).attributes & XMovable) != 0){
                // The cell c is movable and constrained on the right
                auto A = g.addArc(fixed, cell_nodes[c]);
                constraint_arcs.push_back(arc_pair(A, pl.positions_[oc] - pl.widths_[c]));
            }
            else if((circuit.get_cell(oc).attributes & XMovable) != 0){
                // The cell oc is movable and constrained on the left
                auto A = g.addArc(cell_nodes[oc], fixed);
                constraint_arcs.push_back(arc_pair(A, - pl.positions_[c] - pl.widths_[c]));
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
            constraint_arcs.push_back(arc_pair(Ar, pl.max_x_ - pl.widths_[rc]));
        }
    }
    

    // And every pin of every net: arcs too
    for(index_t n=0; n<circuit.net_cnt(); ++n){
        assert(circuit.get_net(n).pin_cnt > 0);
        for(auto p : circuit.get_net(n)){
            index_t c = p.cell_ind;
            int_t pin_offs = 0.5 * pl.widths_[c] + (pl.x_orientations_[c] ? p.offset.x_ : - p.offset.x_); // Offset to the beginning of the cell
            if((circuit.get_cell(c).attributes & XMovable) != 0){
                Arc Al = g.addArc(cell_nodes[c], Lnet_nodes[n]);
                constraint_arcs.push_back(arc_pair(Al, pin_offs));
                Arc Ar = g.addArc(Unet_nodes[n], cell_nodes[c]);
                constraint_arcs.push_back(arc_pair(Ar, -pin_offs));
            }
            else{ // Fixed offset
                auto Al = g.addArc(fixed, Lnet_nodes[n]);
                constraint_arcs.push_back(arc_pair(Al, pl.positions_[c] + pin_offs));
                auto Ar = g.addArc(Unet_nodes[n], fixed);
                constraint_arcs.push_back(arc_pair(Ar, - pl.positions_[c] - pin_offs));
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
    //NetworkSimplex<SmartDigraph> ns(g);
    CostScaling<SmartDigraph, int_t> ns(g);
    //ns.supplyMap(supply).costMap(cost);
    ns.supplyMap(supply).costMap(cost).upperMap(capacity);
    auto res = ns.run();
    std::cout << "Solved the dualMCF problem!" << std::endl;
    if(res == ns.OPTIMAL){
        std::cout << "It is OK!" << std::endl;
    }
    if(res == ns.INFEASIBLE){
        std::cout << "It is unbounded (MCF infeasible)" << std::endl;
        abort();
    }
    if(res == ns.UNBOUNDED){
        std::cout << "It is infeasible (MCF unbounded)" << std::endl;
        abort();
    }
    
    // And we get the new positions as the dual values of the current solution (compared to the fixed pin) 
    for(index_t c=0; c<circuit.cell_cnt(); ++c){ // The cells
        if((circuit.get_cell(c).attributes & XMovable) != 0){
            pl.positions_[c] = ns.potential(cell_nodes[c]) - ns.potential(fixed);
        }
    }
    pl.selfcheck();
}

} // namespace dp
} // namespace coloquinte


