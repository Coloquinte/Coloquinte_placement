#include "Coloquinte/circuit_helper.hxx"

#include <stack>
#include <functional>

namespace coloquinte{
namespace gp{

namespace{
inline void opt_orient(netlist const & circuit, placement_t & pl, std::function<float_t & (point<float_t> &)> coor, mask_t FLIPPABLE){
    std::stack<index_t> opt_cells;
    for(index_t cell_ind = 0; cell_ind < circuit.cell_cnt(); ++cell_ind){
        if( (circuit.get_cell(cell_ind).attributes & FLIPPABLE) != 0) opt_cells.push(cell_ind);
    }
    while(not opt_cells.empty()){
        index_t cell_ind = opt_cells.top(); opt_cells.pop();
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

        if(n_cost < p_cost)
            coor(pl.orientations_[cell_ind]) = -1.0;
        if(p_cost < n_cost)
            coor(pl.orientations_[cell_ind]) =  1.0;

        // If we changed the orientation, check the extreme pins which changed and try their cells again
        if(coor(pl.orientations_[cell_ind]) != old_orientation){
            std::sort(extreme_elements.begin(), extreme_elements.end());
            extreme_elements.resize(std::distance(extreme_elements.begin(), std::unique(extreme_elements.begin(), extreme_elements.end())));
            for(index_t j : extreme_elements){
                opt_cells.push(j);
            }
        }
    }
}
} // End anonymous namespace

// Iteratively optimize feasible orientations; performs only one pass
void optimize_exact_orientations(netlist const & circuit, placement_t & pl){
    opt_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.x_; }, XFlippable);
    opt_orient(circuit, pl, [](point<float_t> & p) -> float_t & { return p.y_; }, YFlippable);
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


