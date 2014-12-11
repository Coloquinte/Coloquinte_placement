/*
 * TODO:
 * Implement something more than basic flipping like
 *      * relaxed optimization (even with LP)
 *      * check if a cell can be fixed immediately, and propagate fixations
 *
 */

#include "Coloquinte/circuit_helper.hxx"

namespace coloquinte{
namespace gp{

// Iteratively optimize feasible orientations; performs only one pass
void optimize_exact_orientations(netlist const & circuit, placement_t & pl){
    for(index_t cell_ind = 0; cell_ind < circuit.cell_cnt(); ++cell_ind){
        if( (circuit.get_cell(cell_ind).attributes & (XFlippable|YFlippable)) == 0) continue; // Not flippable at all, continue

        // Get all the nearby nets
        std::vector<index_t> involved_nets;
        for(netlist::pin_t p : circuit.get_cell(cell_ind)){
            involved_nets.push_back(p.net_ind);
        }
        // Deal with cells with multiple pins in one net (uniquify)
        std::sort(involved_nets.begin(), involved_nets.end());
        involved_nets.resize(std::distance(involved_nets.begin(), std::unique(involved_nets.begin(), involved_nets.end())));

        // Just test every possible orientation and chose the best
        float_t x_p_cost=0.0, x_n_cost=0.0, y_p_cost=0.0, y_n_cost=0.0;
        pl.orientations_[cell_ind] = point<float_t>(1.0, 1.0);
        for(index_t n : involved_nets){
            auto pins = get_pins_1D(circuit, pl, n);
            auto minmaxX = std::minmax_element(pins.x_.begin(), pins.x_.end()),
                 minmaxY = std::minmax_element(pins.y_.begin(), pins.y_.end());
            x_p_cost += minmaxX.second->pos - minmaxX.first->pos;
            y_p_cost += minmaxY.second->pos - minmaxY.first->pos;
        }
        pl.orientations_[cell_ind] = point<float_t>(-1.0, -1.0);
        for(index_t n : involved_nets){
            auto pins = get_pins_1D(circuit, pl, n);
            auto minmaxX = std::minmax_element(pins.x_.begin(), pins.x_.end()),
                 minmaxY = std::minmax_element(pins.y_.begin(), pins.y_.end());
            x_n_cost += minmaxX.second->pos - minmaxX.first->pos;
            y_n_cost += minmaxY.second->pos - minmaxY.first->pos;
        }

        // Test if we have the right to modify this orientation, then chose the best
        if( (circuit.get_cell(cell_ind).attributes & XFlippable) != 0){
            if(x_n_cost < x_p_cost){
                pl.orientations_[cell_ind].x_ = -1.0;
            }
            else{
                pl.orientations_[cell_ind].x_ =  1.0;
            }
        }
        if( (circuit.get_cell(cell_ind).attributes & YFlippable) != 0){
            if(y_n_cost < y_p_cost){
                pl.orientations_[cell_ind].y_ = -1.0;
            }
            else{
                pl.orientations_[cell_ind].y_ =  1.0;
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


