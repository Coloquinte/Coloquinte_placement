
#ifndef COLOQUINTE_GP_CIRCUIT
#define COLOQUINTE_GP_CIRCUIT

#include "common.hxx"
#include "solvers.hxx"
#include "netlist.hxx"
#include "rough_legalizers.hxx"

#include <vector>
#include <cassert>

namespace coloquinte{
namespace gp{

struct placement_t{
    std::vector<point<float_t> > positions_;
    std::vector<point<float_t> > orientations_;

    index_t cell_cnt() const{
        assert(positions_.size() == orientations_.size());
        return positions_.size();
    }
};

point<linear_system> empty_linear_systems(netlist const & circuit, placement_t const & pl);

// Net models stuff
point<linear_system> get_HPWLF_linear_system  (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);
point<linear_system> get_HPWLR_linear_system  (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);
point<linear_system> get_star_linear_system   (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);
point<linear_system> get_clique_linear_system (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);
point<linear_system> get_MST_linear_system    (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);
point<linear_system> get_RSMT_linear_system   (netlist const & circuit, placement_t const & pl, float_t tol, index_t min_s, index_t max_s);

// Additional forces
point<linear_system> get_pulling_forces (netlist const & circuit, placement_t const & pl, float_t typical_distance);
point<linear_system> get_linear_pulling_forces (netlist const & circuit, placement_t const & UB_pl, placement_t const & LB_pl, float_t force, float_t min_distance);

// Solve the final linear system
void get_result(netlist const & circuit, placement_t & pl, point<linear_system> & L, index_t nbr_iter);

// Cost-related stuff, whether wirelength or disruption
float_t get_HPWL_wirelength (netlist const & circuit, placement_t const & pl);
float_t get_MST_wirelength  (netlist const & circuit, placement_t const & pl);
float_t get_RSMT_wirelength (netlist const & circuit, placement_t const & pl);

float_t get_mean_linear_disruption(netlist const & circuit, placement_t const & LB_pl, placement_t const & UB_pl);
float_t get_mean_quadratic_disruption(netlist const & circuit, placement_t const & LB_pl, placement_t const & UB_pl);

// Legalizer-related stuff
region_distribution get_rough_legalizer(netlist const & circuit, placement_t const & pl, box<int_t> surface);
void get_result(netlist const & circuit, placement_t & pl, region_distribution const & legalizer);

// Cell orientation optimization
void optimize_exact_orientations(netlist const & circuit, placement_t & pl);
void zero_orientations(netlist const & circuit, placement_t & pl);
void spread_orientations(netlist const & circuit, placement_t & pl);
void round_orientations(netlist const & circuit, placement_t & pl);

} // namespace gp
} // namespace coloquinte

#endif

