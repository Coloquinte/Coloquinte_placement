
#ifndef COLOQUINTE_GP_HELPERCIRCUIT
#define COLOQUINTE_GP_HELPERCIRCUIT

#include "common.hxx"
#include "circuit.hxx"

namespace coloquinte{
namespace gp{

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
    pin_1D x() const{ return pin_1D(cell_ind, pos.x_, offs.x_, movable); }
    pin_1D y() const{ return pin_1D(cell_ind, pos.y_, offs.y_, movable); }
};

inline float_t dist(pin_2D const a, pin_2D const b){
    point<float_t> diff = a.pos - b.pos;
    return std::abs(diff.x_) + std::abs(diff.y_);
}

inline std::vector<pin_2D>         get_pins_2D(netlist const & circuit, placement_t const & pl, index_t net_ind){
    std::vector<pin_2D> ret;
    for(auto p : circuit.get_net(net_ind)){
        point<float_t> offs = static_cast<point<float_t> >(p.offset) * pl.orientations_[p.cell_ind];
        point<float_t> pos  = static_cast<point<float_t> >(offs)     + pl.positions_[p.cell_ind];

        bool movable = (circuit.get_cell(p.cell_ind).attributes & XMovable) != 0 and (circuit.get_cell(p.cell_ind).attributes & YMovable) != 0;
        ret.push_back(pin_2D(p.cell_ind, pos, offs, movable));
    }
    return ret;
}

inline point<std::vector<pin_1D> > get_pins_1D(netlist const & circuit, placement_t const & pl, index_t net_ind){
    point<std::vector<pin_1D> > ret;
    for(auto p : circuit.get_net(net_ind)){
        point<float_t> offs = static_cast<point<float_t> >(p.offset) * pl.orientations_[p.cell_ind];
        point<float_t> pos  = static_cast<point<float_t> >(offs)     + pl.positions_[p.cell_ind];

        bool x_movable = (circuit.get_cell(p.cell_ind).attributes & XMovable) != 0;
        bool y_movable = (circuit.get_cell(p.cell_ind).attributes & YMovable) != 0;
        ret.x_.push_back(pin_1D(p.cell_ind, pos.x_, offs.x_, x_movable));
        ret.y_.push_back(pin_1D(p.cell_ind, pos.y_, offs.y_, y_movable));
    }
    return ret;
}

} // namespace gp
} // namespace coloquinte

#endif

