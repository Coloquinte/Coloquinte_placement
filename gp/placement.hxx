
#ifndef COLOQUINTE_GP_PLACEMENT
#define COLOQUINTE_GP_PLACEMENT

#include "circuit.hxx"

#include <vector>
#include <memory>


namespace coloquinte{
namespace gp{

struct global_placement{
    struct placement_value{
        float_t x_pos_, y_pos_;
        float_t x_flip_, y_flip_; // Orientation between -1 and +1. Not discrete
    };

    std::shared_ptr<circuit> reference_circuit_;
    std::vector<placement_value> placement_;
};


}; // namespace gp
}; // namespace coloquinte

#endif


