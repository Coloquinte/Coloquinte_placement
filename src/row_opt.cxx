
#include "coloquinte/detailed.hxx"
#include "coloquinte/circuit_helper.hxx"
#include "coloquinte/optimization_subproblems.hxx"

#include <cassert>

namespace coloquinte{
namespace dp{

namespace{
struct net_properties{
    int_t ext_pos_min, ext_pos_max, weight;
    index_t cell_fst, cell_lst;
    std::pair<int_t, int_t> offs_fst, offs_lst;
    bool external_pins;

    net_properties(){
        weight = 1;
        ext_pos_min = std::numeric_limits<int_t>::max();
        ext_pos_max = -std::numeric_limits<int_t>::min();
        cell_fst = std::numeric_limits<index_t>::max();
        cell_lst = 0;
        auto empty_pins = std::pair<int_t, int_t>(ext_pos_min, ext_pos_max);
        offs_fst = empty_pins;
        offs_lst = empty_pins;
        external_pins = false;
    }
};

std::vector<net_properties> get_net_properties(netlist const & circuit, detailed_placement const & pl, std::vector<index_t> const & cells){
    struct order_gettr{
        index_t cell_ind, seq_order;
        bool operator<(order_gettr const o) const{ return cell_ind < o.cell_ind; }
        bool operator<(index_t const o) const{ return cell_ind < o; }
        order_gettr(index_t c, index_t i) : cell_ind(c), seq_order(i) {}
    };

    std::vector<order_gettr> cells_in_row;
    for(index_t i=0; i<cells.size(); ++i){
        cells_in_row.push_back(order_gettr(cells[i],i));
    }
    std::sort(cells_in_row.begin(), cells_in_row.end());

    std::vector<index_t> involved_nets;
    for(index_t c : cells){
        for(netlist::pin_t p : circuit.get_cell(c)){
            involved_nets.push_back(p.net_ind);
        }
    }
    // Uniquify the nets
    std::sort(involved_nets.begin(), involved_nets.end());
    involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

    std::vector<net_properties> net_props;

    for(index_t n : involved_nets){
        net_properties cur_net;
        cur_net.weight = circuit.get_net(n).weight; 

        for(netlist::pin_t p : circuit.get_net(n)){
            auto it = std::lower_bound(cells_in_row.begin(), cells_in_row.end(), p.cell_ind);
            if(it != cells_in_row.end() and it->cell_ind == p.cell_ind){
                // Found a pin on the row
                if(it->seq_order < cur_net.cell_fst){
                    cur_net.cell_fst = it->seq_order;
                    cur_net.offs_fst = std::pair<int_t, int_t>(p.offset.x_, p.offset.x_);
                }
                else if(it->seq_order == cur_net.cell_fst){
                    cur_net.offs_fst.first  = std::min(p.offset.x_, cur_net.offs_fst.first);
                    cur_net.offs_fst.second = std::max(p.offset.x_, cur_net.offs_fst.second);
                }
                if(it->seq_order > cur_net.cell_lst){
                    cur_net.cell_lst = it->seq_order;
                    cur_net.offs_lst = std::pair<int_t, int_t>(p.offset.x_, p.offset.x_);
                }
                else if(it->seq_order == cur_net.cell_lst){
                    cur_net.offs_lst.first  = std::min(p.offset.x_, cur_net.offs_lst.first);
                    cur_net.offs_lst.second = std::max(p.offset.x_, cur_net.offs_lst.second);
                }
            }
            else{ // Found a pin which remains fixed for this round
                cur_net.external_pins = true;

                int_t pos = pl.plt_.positions_[p.cell_ind].x_ + (pl.plt_.orientations_[p.cell_ind].x_ ? p.offset.x_ : circuit.get_cell(p.cell_ind).size.x_ - p.offset.x_);
                cur_net.ext_pos_min = std::min(pos, cur_net.ext_pos_min);
                cur_net.ext_pos_max = std::max(pos, cur_net.ext_pos_max);
            }
        }
        net_props.push_back(cur_net);
    }

    return net_props;
}

// Optimizes an ordered sequence of standard cells on the same row, returns the cost and the corresponding positions
inline std::int64_t optimize_convex_sequence(netlist const & circuit, detailed_placement const & pl, std::vector<index_t> const & cells, std::vector<int_t> & positions, int_t lower_lim, int_t upper_lim){
    struct cell_bound{
        index_t c;
        int_t pos;
        int_t slope;
        bool operator<(cell_bound const o) const{ return c < o.c; }
        cell_bound(index_t order, int_t p, int_t s) : c(order), pos(p), slope(s) {}
    };

    auto net_props = get_net_properties(circuit, pl, cells);

    std::vector<cell_bound> bounds;
    std::vector<int_t> right_slopes(cells.size(), 0);
    for(auto const cur_net : net_props){
        if(cur_net.external_pins){
            index_t fst_c = cells[cur_net.cell_fst],
                    lst_c = cells[cur_net.cell_lst];
            int_t fst_pin_offs = (pl.plt_.orientations_[fst_c].x_ ?
                cur_net.offs_fst.first : circuit.get_cell(fst_c).size.x_ - cur_net.offs_fst.second);
            int_t lst_pin_offs = (pl.plt_.orientations_[lst_c].x_ ?
                cur_net.offs_lst.second : circuit.get_cell(lst_c).size.x_ - cur_net.offs_lst.first);

            bounds.push_back(cell_bound(cur_net.cell_fst, cur_net.ext_pos_min - fst_pin_offs, cur_net.weight));
            bounds.push_back(cell_bound(cur_net.cell_lst, cur_net.ext_pos_max - lst_pin_offs, cur_net.weight));

            right_slopes[cur_net.cell_lst] += cur_net.weight;
        }
        else if(cur_net.cell_fst != cur_net.cell_lst){
            right_slopes[cur_net.cell_lst] += cur_net.weight;
            right_slopes[cur_net.cell_fst] -= cur_net.weight;
        }
    }
    std::sort(bounds.begin(), bounds.end());

    full_single_row OSRP;
    for(index_t i=0, j=0; i<cells.size(); ++i){
        index_t cur_cell_ind = cells[i];

        OSRP.push_cell(circuit.get_cell(cur_cell_ind).size.x_, lower_lim, upper_lim);
        OSRP.push_slope(right_slopes[i]);
        for(; j<bounds.size() and bounds[j].c == i; ++j){
            OSRP.push_bound(bounds[j].pos, bounds[j].slope);
        }
    }

    positions = OSRP.get_placement();

    std::int64_t cost = 0;
    for(auto const cur_net : net_props){
        index_t fst_c = cells[cur_net.cell_fst],
                lst_c = cells[cur_net.cell_lst];
        int_t fst_pin_offs = (pl.plt_.orientations_[fst_c].x_ ?
            cur_net.offs_fst.first : circuit.get_cell(fst_c).size.x_ - cur_net.offs_fst.second);
        int_t lst_pin_offs = (pl.plt_.orientations_[lst_c].x_ ?
            cur_net.offs_lst.second : circuit.get_cell(lst_c).size.x_ - cur_net.offs_lst.first);

        if(cur_net.external_pins){
            cost += (std::max(cur_net.ext_pos_max, positions[cur_net.cell_lst] + lst_pin_offs) - std::min(cur_net.ext_pos_min, positions[cur_net.cell_fst] + fst_pin_offs));
        }
        else{
            cost += (lst_pin_offs - fst_pin_offs);
        }
    }
    return cost;
}
} // End anonymous namespace

void OSRP_convex_HPWL(netlist const & circuit, detailed_placement & pl){
    for(index_t r=0; r<pl.row_cnt(); ++r){
        index_t OSRP_cell = pl.get_first_cell_on_row(r);

        while(OSRP_cell != null_ind){
            std::vector<index_t> cells;

            for(;
                OSRP_cell != null_ind and pl.cell_height(OSRP_cell) == 1 and (circuit.get_cell(OSRP_cell).attributes & XMovable) != 0; // Movable standard cell
                OSRP_cell = pl.neighbours_[pl.neighbours_limits_[OSRP_cell]].second
            ){
                cells.push_back(OSRP_cell);
            }

            if(not cells.empty()){
                // Limits of the placement region for the cells taken
                int_t lower_lim = pl.get_limit_positions(circuit, cells.front()).first,
                      upper_lim = pl.get_limit_positions(circuit, cells.back()).second;
                std::vector<int_t> final_positions;
                optimize_convex_sequence(circuit, pl, cells, final_positions, lower_lim, upper_lim);

                // Update the positions
                for(index_t i=0; i<cells.size(); ++i){
                    pl.plt_.positions_[cells[i]].x_ = final_positions[i];
                }
            }

            if(OSRP_cell != null_ind) OSRP_cell = pl.get_next_cell_on_row(OSRP_cell, r); // Go to the next group
        } // Iteration on the entire row
    } // Iteration on the rows

    pl.selfcheck();
}

void swaps_row_HPWL(netlist const & circuit, detailed_placement & pl, index_t range){
    assert(range >= 2);

    for(index_t r=0; r<pl.row_cnt(); ++r){
        index_t OSRP_cell = pl.get_first_cell_on_row(r);

        while(OSRP_cell != null_ind){
            std::vector<index_t> cells;
            for(index_t nbr_cells=0;
                    OSRP_cell != null_ind
                and pl.cell_height(OSRP_cell) == 1
                and (circuit.get_cell(OSRP_cell).attributes & XMovable) != 0
                and nbr_cells < range;
                OSRP_cell = pl.neighbours_[pl.neighbours_limits_[OSRP_cell]].second, ++nbr_cells
            ){
                cells.push_back(OSRP_cell);
            }

            if(not cells.empty()){
                // Limits of the placement region for the cells taken
                int_t lower_lim = pl.get_limit_positions(circuit, cells.front()).first,
                      upper_lim = pl.get_limit_positions(circuit, cells.back()).second;

                std::int64_t best_cost = std::numeric_limits<std::int64_t>::max();
                std::vector<int_t> positions(cells.size());
                std::vector<index_t> best_permutation = cells;
                std::vector<int_t> best_positions(cells.size());;
                std::vector<index_t> old_permutation = cells;

                // Check every possible permutation of the cells
                std::sort(cells.begin(), cells.end());
                do{
                    std::int64_t cur_cost = optimize_convex_sequence(circuit, pl, cells, positions, lower_lim, upper_lim);
                    assert(std::isfinite(cur_cost) );
                    if(cur_cost < best_cost){
                        best_cost = cur_cost;
                        best_permutation = cells;
                        best_positions = positions;
                    }
                }while(std::next_permutation(cells.begin(), cells.end()));
    
                cells = best_permutation;
                // Update the positions and the topology
                for(index_t i=0; i<cells.size(); ++i){
                    pl.plt_.positions_[cells[i]].x_ = best_positions[i];
                }

                pl.reorder_standard_cells(old_permutation, cells);
           }
    
            if(OSRP_cell != null_ind){
                // We are on a non-movable cell
                if( (circuit.get_cell(OSRP_cell).attributes & XMovable) == 0 or pl.cell_height(OSRP_cell) != 1){
                    OSRP_cell = pl.get_next_cell_on_row(OSRP_cell, r); // Go to the next group
                }
                else{ // We optimized with the maximum number of cells: just advance a few cells and optimize again
                    assert(cells.size() == range);
                    OSRP_cell = cells[range/2];
                }
            }
        } // Iteration on the entire row
    } // Iteration on the rows

    pl.selfcheck();
}

} // namespace dp
} // namespace coloquinte


