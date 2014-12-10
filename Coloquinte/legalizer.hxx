
#include "circuit.hxx"
#include "detailed.hxx"

namespace coloquinte{
namespace dp{

detailed_placement legalize(netlist const & circuit, gp::placement_t const & pl, box<int_t> surface, int_t row_height);

} // namespace dp
} // namespace coloquinte
