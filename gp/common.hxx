
#ifndef COLOQUINTE_GP_COMMON
#define COLOQUINTE_GP_COMMON

#include <cstdint>

namespace coloquinte{

using float_t    = float;
using int_t      = std::int32_t;
using index_t    = std::uint32_t;
using capacity_t = std::int64_t;
using mask_t     = std::uint32_t;

using ext_object = std::uint64_t;

enum PlacementType{
    Optimist  = 0,
    Pessimist = 1
};

enum Movability{
    XMovable   = 1     ,
    YMovable   = 1 << 1,
    XFlippable = 1 << 2,
    YFlippable = 1 << 3,
    SoftMacro  = 1 << 4
};

template<typename T>
struct point{
    T x_, y_;
    point(){}
    point(T x, T y): x_(x), y_(y){}

    template<typename S>
    point<S> cast() const{
        return point<S>(static_cast<S>(x_), static_cast<S>(y_));
    }
    template<typename S>
    operator point<S>(){
        return cast<S>();
    }

    point<T> operator+(point<T> const o) const{
        return point<T>(x_+o.x_, y_+o.y_);
    }
    point<T> operator-(point<T> const o) const{
        return point<T>(x_-o.x_, y_-o.y_);
    }
    point<T> scale(T const o) const{
        return point<T>(o*x_, o*y_);
    }
};

template<typename T>
struct box{
    T x_min_, y_min_, x_max_, y_max_;
    box(){}
    box(T x_mn, T x_mx, T y_mn, T y_mx) : x_min_(x_mn), x_max_(x_mx), y_min_(y_mn), y_max_(y_mx){}
};

using orientation_t = point<bool>;

} // Namespace coloquinte

#endif

