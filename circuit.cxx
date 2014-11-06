#include "gp/circuit.hxx"

#include <cassert>

namespace coloquinte{
namespace gp{

void add_forces(circuit::pin_1D const p1, circuit::pin_1D const p2, linear_system & L, float_t tol, float_t scale){
    if(p1.movable && p2.movable){
        L.add_force(
            tol,
            scale,
            p1.cell_ind, p2.cell_ind,
            p1.pos,      p2.pos,
            p1.offs,     p2.offs
        );
    }
    else if(p1.movable){
        L.add_fixed_force(
            tol, scale,
            p1.cell_ind,
            p1.pos, p2.pos,
            p1.offs
        );
    }
    else if(p2.movable){
        L.add_fixed_force(
            tol, scale,
            p2.cell_ind,
            p2.pos, p1.pos,
            p2.offs
        );
    }
}

circuit::circuit(box<int_t> placement_surface, std::vector<point<float_t> > cell_positions, std::vector<point<float_t> > cell_orientations,
                 std::vector<temporary_cell> all_cells, std::vector<temporary_net> all_nets, std::vector<temporary_pin> all_pins)
            : internal_netlist(all_cells, all_nets, all_pins), placement_area_(placement_surface){
    placement_t common;
    common.positions_ = cell_positions;
    common.orientations_ = cell_orientations;
    placements_.resize(2, common);
}

point<linear_system> circuit::empty_linear_systems() const{
    point<linear_system> ret = point<linear_system>(linear_system(cell_cnt()), linear_system(cell_cnt()));

    /*
    for(index_t i=0; i<cell_cnt(); ++i){
        if( (XMovable & internal_netlist.get_cell(i).attributes) == 0){
            ret.x_.add_triplet(i, i, 1.0);
        }
        if( (YMovable & internal_netlist.get_cell(i).attributes) == 0){
            ret.y_.add_triplet(i, i, 1.0);
        }
    }
    */

    return ret;
}

std::vector<circuit::pin_2D>         circuit::get_pins_2D(index_t placement_ind, index_t net_ind) const{
    std::vector<pin_2D> ret;
    for(auto p : internal_netlist.get_net(net_ind)){
        point<float_t> offs = static_cast<point<float_t> >(p.offset) * placements_[placement_ind].orientations_[p.cell_ind];
        point<float_t> pos  = static_cast<point<float_t> >(offs)     + placements_[placement_ind].positions_[p.cell_ind];

        bool movable = (internal_netlist.get_cell(p.cell_ind).attributes & (XMovable|YMovable)) != 0;
        ret.push_back(pin_2D(p.cell_ind, pos, offs, movable));
    }
    return ret;
}
point<std::vector<circuit::pin_1D> > circuit::get_pins_1D(index_t placement_ind, index_t net_ind) const{
    point<std::vector<circuit::pin_1D> > ret;
    for(auto p : internal_netlist.get_net(net_ind)){
        point<float_t> offs = static_cast<point<float_t> >(p.offset) * placements_[placement_ind].orientations_[p.cell_ind];
        point<float_t> pos  = static_cast<point<float_t> >(offs)     + placements_[placement_ind].positions_[p.cell_ind];

        bool x_movable = (internal_netlist.get_cell(p.cell_ind).attributes & XMovable) != 0;
        bool y_movable = (internal_netlist.get_cell(p.cell_ind).attributes & YMovable) != 0;
        ret.x_.push_back(pin_1D(p.cell_ind, pos.x_, offs.x_, x_movable));
        ret.y_.push_back(pin_1D(p.cell_ind, pos.y_, offs.y_, y_movable));
    }
    return ret;
}

void get_HPWLF(std::vector<circuit::pin_1D> const & pins, linear_system & L, float_t tol){
    if(pins.size() >= 2){
        auto min_elt = std::min_element(pins.begin(), pins.end()), max_elt = std::max_element(pins.begin(), pins.end());

        for(auto it = pins.begin(); it != pins.end(); ++it){
            if(it != min_elt){
                add_forces(*it, *min_elt, L, tol, 1.0/(pins.size()-1));
            }
            else if(it != max_elt){
                add_forces(*it, *max_elt, L, tol, 1.0/(pins.size()-1));
            }
        }
    }
}

float_t circuit::get_HPWL_wirelength(index_t placement_ind) const{
    float_t sum = 0.0;
    for(index_t i=0; i<net_cnt(); ++i){
        auto pins = get_pins_1D(placement_ind, i);
        auto minmaxX = std::minmax_element(pins.x_.begin(), pins.x_.end()), minmaxY = std::minmax_element(pins.y_.begin(), pins.y_.end());
        sum += ((minmaxX.second->pos - minmaxX.first->pos) + (minmaxY.second->pos - minmaxY.first->pos));
    }
    return sum;
}

point<linear_system> circuit::get_HPWLF_linear_system (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const{
    point<linear_system> L = empty_linear_systems();
    for(index_t i=0; i<net_cnt(); ++i){
        auto pins = get_pins_1D(placement_ind, i);
        get_HPWLF(pins.x_, L.y_, tol);
        get_HPWLF(pins.y_, L.y_, tol);
    }
    return L;
}

void get_HPWLR(std::vector<circuit::pin_1D> const & pins, linear_system & L, float_t tol){
    std::vector<circuit::pin_1D> sorted_pins = pins;
    std::sort(sorted_pins.begin(), sorted_pins.end());
    for(index_t i=0; i+1<sorted_pins.size(); ++i){
        add_forces(sorted_pins[i], sorted_pins[i+1], L, tol, 1.0);
    }
}

point<linear_system> circuit::get_HPWLR_linear_system (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const{
    point<linear_system> L = empty_linear_systems();
    for(index_t i=0; i<net_cnt(); ++i){
        auto pins = get_pins_1D(placement_ind, i);
        get_HPWLR(pins.x_, L.y_, tol);
        get_HPWLR(pins.y_, L.y_, tol);
    }
    return L;
}


void circuit::get_result(float_t tol, index_t placement_ind, point<linear_system> & L){
    std::vector<float_t> x_sol, y_sol;
    std::vector<float_t> x_guess(cell_cnt()), y_guess(cell_cnt());
    
    assert(L.x_.size() == x_guess.size());
    assert(L.y_.size() == y_guess.size());

    for(index_t i=0; i<cell_cnt(); ++i){
        x_guess[i] = placements_[placement_ind].positions_[i].x_;
        y_guess[i] = placements_[placement_ind].positions_[i].y_;
    }
    
    #pragma omp task shared(L) shared(x_sol)
    x_sol = L.x_.solve_CG(x_guess, tol);
    //x_sol = L.x_.solve_cholesky();
    #pragma omp task shared(L) shared(y_sol)
    y_sol = L.y_.solve_CG(y_guess, tol);
    //y_sol = L.y_.solve_cholesky();
    #pragma omp taskwait
    
    for(index_t i=0; i<cell_cnt(); ++i){
        // TODO: correct formulation
        if( (internal_netlist.get_cell(i).attributes & (XMovable|YMovable)) != 0){
            placements_[placement_ind].positions_[i] = point<float_t>(x_sol[i], y_sol[i]);
        }
    }
}

region_distribution circuit::get_rough_legalizer(index_t placement_ind) const{
    std::vector<region_distribution::movable_cell> movable_cells;
    std::vector<region_distribution::fixed_cell>   fixed_cells;

    for(index_t i=0; i<cell_cnt(); ++i){
        auto C = internal_netlist.get_cell(i);
        if(C.attributes & (XMovable|YMovable)){
            movable_cells.push_back(region_distribution::movable_cell(C.area, placements_[placement_ind].positions_[i], i));
        }
        else{
            fixed_cells.push_back(region_distribution::fixed_cell(C.size, placements_[placement_ind].positions_[i]));
        }
    }

    return region_distribution(placement_area_, movable_cells, fixed_cells);
}


} // namespace gp
} // namespace coloquinte


