#include "coloquinte/circuit_helper.hxx"

#include <stack>
#include <functional>
#include <algorithm>

namespace coloquinte{
namespace gp{

namespace{
index_t const null_ind = std::numeric_limits<index_t>::max();
float_t const INF = std::numeric_limits<float_t>::infinity();

inline void opt_orient(netlist const & circuit, placement_t & pl, std::function<float_t & (point<float_t> &)> coor, mask_t FLIPPABLE){
    std::stack<index_t> opt_cells;
    for(index_t cell_ind = 0; cell_ind < circuit.cell_cnt(); ++cell_ind){
        if( (circuit.get_cell(cell_ind).attributes & FLIPPABLE) != 0)
            opt_cells.push(cell_ind);
    }
    while(not opt_cells.empty()){
        index_t cell_ind = opt_cells.top(); opt_cells.pop();
        assert((circuit.get_cell(cell_ind).attributes & FLIPPABLE) != 0);

        // What is the current orientation?
        float_t old_orientation = coor(pl.orientations_[cell_ind]);
        float_t pos = coor(pl.positions_[cell_ind]);

        // Check both orientations of the cell
        std::vector<index_t> involved_nets;
        for(netlist::pin_t p : circuit.get_cell(cell_ind)){
            involved_nets.push_back(p.net_ind);
        }
        // Deal with cells with multiple pins in one net (uniquify)
        std::sort(involved_nets.begin(), involved_nets.end());
        involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

        float_t p_cost = 0.0, n_cost = 0.0;
        std::vector<index_t> extreme_elements;
        for(index_t n : involved_nets){
            std::vector<pin_1D> other_pins;
            std::vector<float_t> offsets;
            for(auto p : circuit.get_net(n)){
                if(p.cell_ind != cell_ind){
                    other_pins.push_back(pin_1D(
                        p.cell_ind,
                        coor(pl.positions_[p.cell_ind]) + coor(p.offset) * coor(pl.orientations_[p.cell_ind]),
                        std::numeric_limits<float_t>::quiet_NaN(), // Don't care about the offset
                        (circuit.get_cell(p.cell_ind).attributes & FLIPPABLE) != 0)
                    );
                }
                else{
                    offsets.push_back(coor(p.offset));
                }
            }
            assert(offsets.size() > 0);
            if(other_pins.size() > 0){ // Else the orientation of the cell doesn't change anything
                auto minmaxC = std::minmax_element(other_pins.begin(), other_pins.end());
                auto minmaxO = std::minmax_element(offsets.begin(), offsets.end());
                p_cost += std::max(pos + *minmaxO.second, minmaxC.second->pos) - std::min(pos + *minmaxO.first, minmaxC.first->pos);
                n_cost += std::max(pos - *minmaxO.first, minmaxC.second->pos) - std::min(pos - *minmaxO.second, minmaxC.first->pos);

                // Do the extreme elements change between the two positions?
                if(minmaxC.second->movable
               and minmaxC.second->pos < std::max(pos + *minmaxO.second, pos - *minmaxO.first)
               and minmaxC.second->pos > std::min(pos + *minmaxO.second, pos - *minmaxO.first)){
                    extreme_elements.push_back(minmaxC.second->cell_ind);
                }
                if(minmaxC.first->movable
               and minmaxC.first->pos < std::max(pos - *minmaxO.second, pos + *minmaxO.first)
               and minmaxC.first->pos > std::min(pos - *minmaxO.second, pos + *minmaxO.first)){
                    extreme_elements.push_back(minmaxC.first->cell_ind);
                }
            }
        }

        coor(pl.orientations_[cell_ind]) = p_cost <= n_cost ? 1.0 : -1.0;

        // If we changed the orientation, check the extreme pins which changed and try their cells again
        if(coor(pl.orientations_[cell_ind]) != old_orientation){
            std::sort(extreme_elements.begin(), extreme_elements.end());
            extreme_elements.resize(std::distance(extreme_elements.begin(), std::unique(extreme_elements.begin(), extreme_elements.end())));
            for(index_t extreme_cell : extreme_elements){
                if( (circuit.get_cell(extreme_cell).attributes & FLIPPABLE) != 0)
                    opt_cells.push(extreme_cell);
            }
        }
    }
}

inline void spread_orient(netlist const & circuit, placement_t & pl, std::function<float_t & (point<float_t> &)> coor, mask_t FLIPPABLE){
    std::vector<float_t> weights(circuit.cell_cnt(), 0.0);
    for(index_t n=0; n<circuit.net_cnt(); ++n){
        float_t min_pos=INF, max_pos=-INF;
        float_t min_offs=INF, max_offs=-INF;
        index_t min_ind=null_ind, max_ind=null_ind;
        for(netlist::pin_t p : circuit.get_net(n)){
            if( (circuit.get_cell(p.cell_ind).attributes & FLIPPABLE) != 0){
                float_t pos = coor(pl.positions_[p.cell_ind]);
                if(pos < min_pos){
                    min_pos = pos;
                    min_ind = p.cell_ind;
                    min_offs = coor(p.offset);
                }
                if(pos > max_pos){
                    max_pos = pos;
                    max_ind = p.cell_ind;
                    max_offs = coor(p.offset);
                }
            }
            else{
                float_t pos = coor(pl.positions_[p.cell_ind]) + coor(pl.orientations_[p.cell_ind]) * coor(p.offset);
                if(pos < min_pos){
                    min_pos = pos;
                    min_ind = null_ind;
                }
                if(pos > max_pos){
                    max_pos = pos;
                    max_ind = null_ind;
                }
            }
        }

        float_t net_weight = circuit.get_net(n).weight;

        if(min_ind != null_ind) weights[min_ind] += net_weight * min_offs;
        if(max_ind != null_ind) weights[max_ind] -= net_weight * max_offs;
    }

    for(index_t c=0; c<circuit.cell_cnt(); ++c){
        coor(pl.orientations_[c]) = (weights[c] >= 0.0) ? 1.0 : -1.0;
    }
}

} // End anonymous namespace

// Iteratively optimize feasible orientations; performs only one pass
void optimize_exact_orientations(netlist const & circuit, placement_t & pl){
    opt_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.x_; }, XFlippable);
    opt_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.y_; }, YFlippable);
}

void spread_orientations(netlist const & circuit, placement_t & pl){
    spread_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.x_; }, XFlippable);
    spread_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.y_; }, YFlippable);
}

void zero_orientations(netlist const & circuit, placement_t & pl){
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        if( (circuit.get_cell(i).attributes & XFlippable) != 0){
            pl.orientations_[i].x_ =  0.0;
        }
        if( (circuit.get_cell(i).attributes & YFlippable) != 0){
            pl.orientations_[i].y_ =  0.0;
        }
    }
}



// Just get the feasible orientation closest to the current solution
void round_orientations(netlist const & circuit, placement_t & pl){
    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        if( (circuit.get_cell(i).attributes & XFlippable) != 0){
            if(pl.orientations_[i].x_ >= 0.0){
                pl.orientations_[i].x_ =  1.0;
            }
            else{
                pl.orientations_[i].x_ = -1.0;
            }
        }

        if( (circuit.get_cell(i).attributes & YFlippable) != 0){
            if(pl.orientations_[i].y_ >= 0.0){
                pl.orientations_[i].y_ =  1.0;
            }
            else{
                pl.orientations_[i].y_ = -1.0;
            }
        }
    }
}

} // namespace gp
} // namespace coloquinte


