/*
 * TODO:
 * Implement something more than basic flipping like
 *      * relaxed optimization (even with LP)
 *      * check if a cell can be fixed immediately, and propagate fixations
 *
 */

#include "Coloquinte/circuit_helper.hxx"

#include <stack>

namespace coloquinte{
namespace gp{

// Iteratively optimize feasible orientations; performs only one pass
void optimize_exact_orientations(netlist const & circuit, placement_t & pl){
    std::stack<index_t> x_opt_cells, y_opt_cells;
    for(index_t cell_ind = 0; cell_ind < circuit.cell_cnt(); ++cell_ind){
        if( (circuit.get_cell(cell_ind).attributes & XFlippable) != 0) x_opt_cells.push(cell_ind);
        if( (circuit.get_cell(cell_ind).attributes & YFlippable) != 0) y_opt_cells.push(cell_ind);
    }

    while(not x_opt_cells.empty()){
        index_t cell_ind = x_opt_cells.top(); x_opt_cells.pop();
        // What is the current orientation?
        float_t old_orientation = pl.orientations_[cell_ind].x_;
        float_t pos = pl.positions_[cell_ind].x_;

        // Check both orientations of the cell
        std::vector<index_t> involved_nets;
        for(netlist::pin_t p : circuit.get_cell(cell_ind)){
            involved_nets.push_back(p.net_ind);
        }
        // Deal with cells with multiple pins in one net (uniquify)
        std::sort(involved_nets.begin(), involved_nets.end());
        involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

        float_t x_p_cost = 0.0, x_n_cost = 0.0;
        std::vector<index_t> extreme_elements;
        for(index_t n : involved_nets){
            std::vector<pin_1D> other_pins;
            std::vector<float_t> offsets;
            for(auto p : circuit.get_net(n)){
                if(p.cell_ind != cell_ind){
                    other_pins.push_back(pin_1D(
                        p.cell_ind,
                        pl.positions_[p.cell_ind].x_ + p.offset.x_ * pl.orientations_[p.cell_ind].x_,
                        std::numeric_limits<float_t>::quiet_NaN(), // Don't care about the offset
                        (circuit.get_cell(p.cell_ind).attributes & XMovable) != 0)
                    );
                }
                else{
                    offsets.push_back(p.offset.x_);
                }
            }
            assert(offsets.size() > 0);
            if(other_pins.size() > 0){ // Else the orientation of the cell doesn't change anything
                auto minmaxX = std::minmax_element(other_pins.begin(), other_pins.end());
                auto minmaxO = std::minmax_element(offsets.begin(), offsets.end());
                x_p_cost += std::max(pos + *minmaxO.second, minmaxX.second->pos) - std::min(pos + *minmaxO.first, minmaxX.first->pos);
                x_n_cost += std::max(pos - *minmaxO.first, minmaxX.second->pos) - std::min(pos - *minmaxO.second, minmaxX.first->pos);

                // Do the extreme elements change between the two positions?
                if(minmaxX.second->movable
               and minmaxX.second->pos < std::max(pos + *minmaxO.second, pos - *minmaxO.first)
               and minmaxX.second->pos > std::min(pos + *minmaxO.second, pos - *minmaxO.first)){
                    extreme_elements.push_back(minmaxX.second->cell_ind);
                }
                if(minmaxX.first->movable
               and minmaxX.first->pos < std::max(pos - *minmaxO.second, pos + *minmaxO.first)
               and minmaxX.first->pos > std::min(pos - *minmaxO.second, pos + *minmaxO.first)){
                    extreme_elements.push_back(minmaxX.first->cell_ind);
                }
            }
        }

        if(x_n_cost < x_p_cost)
            pl.orientations_[cell_ind].x_ = -1.0;
        if(x_p_cost < x_n_cost)
            pl.orientations_[cell_ind].x_ =  1.0;

        // If we changed the orientation, check the extreme pins which changed and try their cells again
        if(pl.orientations_[cell_ind].x_ != old_orientation){
            std::sort(extreme_elements.begin(), extreme_elements.end());
            extreme_elements.resize(std::distance(extreme_elements.begin(), std::unique(extreme_elements.begin(), extreme_elements.end())));
            for(index_t j : extreme_elements){
                x_opt_cells.push(j);
            }
        }
    }

    while(not y_opt_cells.empty()){
        index_t cell_ind = y_opt_cells.top(); y_opt_cells.pop();
        // What is the current orientation?
        float_t old_orientation = pl.orientations_[cell_ind].y_;
        float_t pos = pl.positions_[cell_ind].y_;

        // Check both orientations of the cell
        std::vector<index_t> involved_nets;
        for(netlist::pin_t p : circuit.get_cell(cell_ind)){
            involved_nets.push_back(p.net_ind);
        }
        // Deal with cells with multiple pins in one net (uniquify)
        std::sort(involved_nets.begin(), involved_nets.end());
        involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

        float_t y_p_cost = 0.0, y_n_cost = 0.0;
        std::vector<index_t> extreme_elements;
        for(index_t n : involved_nets){
            std::vector<pin_1D> other_pins;
            std::vector<float_t> offsets;
            for(auto p : circuit.get_net(n)){
                if(p.cell_ind != cell_ind){
                    other_pins.push_back(pin_1D(
                        p.cell_ind,
                        pl.positions_[p.cell_ind].y_ + p.offset.y_ * pl.orientations_[p.cell_ind].y_,
                        std::numeric_limits<float_t>::quiet_NaN(), // Don't care about the offset
                        (circuit.get_cell(p.cell_ind).attributes & YMovable) != 0)
                    );
                }
                else{
                    offsets.push_back(p.offset.y_);
                }
            }
            assert(offsets.size() > 0);
            if(other_pins.size() > 0){ // Else the orientation of the cell doesn't change anything
                auto minmaxY = std::minmax_element(other_pins.begin(), other_pins.end());
                auto minmaxO = std::minmax_element(offsets.begin(), offsets.end());
                y_p_cost += std::max(pos + *minmaxO.second, minmaxY.second->pos) - std::min(pos + *minmaxO.first, minmaxY.first->pos);
                y_n_cost += std::max(pos - *minmaxO.first, minmaxY.second->pos) - std::min(pos - *minmaxO.second, minmaxY.first->pos);

                // Do the extreme elements change between the two positions?
                if(minmaxY.second->movable
               and minmaxY.second->pos < std::max(pos + *minmaxO.second, pos - *minmaxO.first)
               and minmaxY.second->pos > std::min(pos + *minmaxO.second, pos - *minmaxO.first)){
                    extreme_elements.push_back(minmaxY.second->cell_ind);
                }
                if(minmaxY.first->movable
               and minmaxY.first->pos < std::max(pos - *minmaxO.second, pos + *minmaxO.first)
               and minmaxY.first->pos > std::min(pos - *minmaxO.second, pos + *minmaxO.first)){
                    extreme_elements.push_back(minmaxY.first->cell_ind);
                }
            }
        }

        if(y_n_cost < y_p_cost)
            pl.orientations_[cell_ind].y_ = -1.0;
        if(y_p_cost < y_n_cost)
            pl.orientations_[cell_ind].y_ =  1.0;

        // If we changed the orientation, check the extreme pins which changed and try their cells again
        if(pl.orientations_[cell_ind].y_ != old_orientation){
            std::sort(extreme_elements.begin(), extreme_elements.end());
            extreme_elements.resize(std::distance(extreme_elements.begin(), std::unique(extreme_elements.begin(), extreme_elements.end())));
            for(index_t j : extreme_elements){
                y_opt_cells.push(j);
            }
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


