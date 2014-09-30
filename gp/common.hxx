
#ifndef COLOQUINTE_GP_COMMON
#define COLOQUINTE_GP_COMMON

#include <cstdint>

namespace coloquinte{
namespace gp{

using float_t = float;
using int_t   = std::int32_t;
using index_t = std::uint32_t;
using capacity_t = std::int64_t;

using ext_object = std::uint64_t;

template<typename T>
struct box{
    T x_min_, y_min_, x_max_, y_max_;
    box(){}
    box(T x_mn, T x_mx, T y_mn, T y_mx) : x_min_(x_mn), x_max_(x_mx), y_min_(y_mn), y_max_(y_mx){}
};

} // Namespace gp
} // Namespace coloquinte

#endif

