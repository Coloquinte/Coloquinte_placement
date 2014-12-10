
#include "circuit.hxx"
#include "detailed.hxx"

namespace coloquinte{
namespace dp{

detailed_placement legalize(netlist const & circuit, gp::placement_t const & pl, box<int_t> surface, int_t row_height);

float_t get_mean_linear_disruption(netlist const & circuit, gp::placement_t const & gpl, detailed_placement const & dpl);
float_t get_mean_quadratic_disruption(netlist const & circuit, gp::placement_t const & gpl, detailed_placement const & dpl);

} // namespace dp
} // namespace coloquinte
