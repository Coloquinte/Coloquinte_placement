
#ifndef COLOQUINTE_GP_CIRCUIT
#define COLOQUINTE_GP_CIRCUIT

#include "common.hxx"
#include "solvers.hxx"
#include "netlist.hxx"
#include "rough_legalizers.hxx"

#include <vector>


namespace coloquinte{

namespace gp{

// Main class
struct circuit{
    struct placement_t{
        std::vector<point<float_t> > positions_;
        std::vector<point<float_t> > orientations_;
    };

    // Members
    coloquinte::netlist internal_netlist;

    box<int_t>                 placement_area_;
    std::vector<placement_t>   placements_; // Current placements (at least two: one optimistic at index 0 and one pessimistic at index 1)

    // Helpers
    public:

    struct pin_1D{
        index_t cell_ind;
        float_t pos;
        float_t offs;
        bool movable;

        bool operator<(pin_1D const o) const { return pos < o.pos; }

        pin_1D(index_t c, float_t p, float_t o, bool m) : cell_ind(c), pos(p), offs(o), movable(m){}
    };
    struct pin_2D{
        index_t        cell_ind;
        point<float_t> pos;
        point<float_t> offs;
        bool movable;

        pin_2D(index_t c, point<float_t> p, point<float_t> o, bool m) : cell_ind(c), pos(p), offs(o), movable(m){}
    };

    

    private:
    point<std::vector<pin_1D> > get_pins_1D(index_t placement_ind, index_t net_ind) const;
    std::vector<pin_2D>         get_pins_2D(index_t placement_ind, index_t net_ind) const;

    point<linear_system> empty_linear_systems() const;

    public:
    circuit(box<int_t> placement_surface, std::vector<point<float_t> > cell_positions, std::vector<point<float_t> > cell_orientations, std::vector<temporary_cell> all_cells, std::vector<temporary_net> all_nets, std::vector<temporary_pin> all_pins);


    point<linear_system> get_HPWLF_linear_system (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const;
    point<linear_system> get_HPWLR_linear_system (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const;
    point<linear_system> get_star_linear_system  (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const;
    point<linear_system> get_MST_linear_system   (float_t tol, index_t placement_ind, index_t min_s, index_t max_s) const;

    point<linear_system> get_pulling_forces (float_t strength, index_t placement_ind) const;

    float_t get_HPWL_wirelength(index_t placement_ind) const;

    region_distribution get_rough_legalizer() const;

    index_t cell_cnt() const{ return internal_netlist.cell_cnt(); }
    index_t net_cnt() const{ return internal_netlist.net_cnt(); }
};



} // namespace gp
} // namespace coloquinte

#endif

