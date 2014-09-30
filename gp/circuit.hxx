
#ifndef COLOQUINTE_GP_MODULES
#define COLOQUINTE_GP_MODULES

#include "common.hxx"
#include <vector>

namespace coloquinte{
namespace gp{

struct circuit{
    public:
    // Nets, their weights and their pins
    struct net_pin{
        index_t cell_ind_, pin_ind_in_cell_;
    };
    struct net{
        index_t nbr_pins_;
        float_t weight_;
        net_pin * pins_;
        ext_object backpointer_;

        net(index_t pin_cnt);
        ~net();
    };

    // Cells, their areas, their pins and their mobility attributes
    struct cell_pin{
        index_t net_ind_, pin_ind_in_net_;
        index_t x_offset_, y_offset_;
    };
    struct cell{
        class attributes{
            // Movability values
            static const std::uint32_t Movable = 1;
            static const std::uint32_t XFlippable = 1 << 1;
            static const std::uint32_t YFlippable = 1 << 2;
            static const std::uint32_t SoftMacro = 1 << 3;
        
            std::uint32_t attr_;
        
            public:
            attributes();

            bool is_fully_movable() const{ return (attr_ & Movable) != 0; }
            bool is_x_flippable()   const{ return (attr_ & XFlippable) != 0; }
            bool is_y_flippable()   const{ return (attr_ & YFlippable) != 0; }
        
            void set_fully_movable(bool mv){  attr_ = mv ? attr_ | Movable    : attr_ & ~Movable;}
            void set_x_flippable  (bool mv){  attr_ = mv ? attr_ | XFlippable : attr_ & ~XFlippable;}
            void set_y_flippable  (bool mv){  attr_ = mv ? attr_ | YFlippable : attr_ & ~YFlippable;}
        };

        capacity_t area_;
        float_t x_size_, y_size_;

        attributes movability_;

        index_t nbr_pins_;
        cell_pin * pins_;
        ext_object backpointer_;

        cell(index_t pin_cnt);
        ~cell();
    };

    private:
    // Members
    box<int_t> placement_area_;
    std::vector<cell> cells_;
    std::vector<net> nets_;


    public:
    inline index_t cell_cnt() const;
    inline index_t net_cnt() const;

    // Getters used by most modules (not meant to modify the circuit)
    inline cell const & get_cell(index_t index) const;
    inline net  const & get_net(index_t index)  const;

    inline box<int_t> placement_area() const;

    // Net modifiers
    inline void set_weight(index_t index, float_t weight);

    // Cell modifiers
    inline void set_area(index_t index, capacity_t area);
    inline void set_size(index_t index, float_t x_size, float_t y_size);
};

inline circuit::net::net(index_t pin_cnt) : nbr_pins_(pin_cnt), pins_(new net_pin[pin_cnt]){}
inline circuit::net::~net(){delete pins_;}

inline circuit::cell::attributes::attributes() : attr_(0){}

inline circuit::cell::cell(index_t pin_cnt) : nbr_pins_(pin_cnt), pins_(new cell_pin[pin_cnt]){}
inline circuit::cell::~cell(){delete pins_;}

inline index_t circuit::cell_cnt() const{ return cells_.size(); }
inline index_t circuit::net_cnt() const{ return nets_.size(); }
inline circuit::cell const & circuit::get_cell(index_t index) const{ return cells_[index]; }
inline circuit::net  const & circuit::get_net(index_t index)  const{ return nets_[index]; }

inline void circuit::set_weight(index_t index, float_t weight){
    nets_[index].weight_ = weight; 
}

inline box<int_t> circuit::placement_area() const{
    return placement_area_;
}

inline void circuit::set_area(index_t index, capacity_t area){
    cells_[index].area_ = area;
}
inline void circuit::set_size(index_t index, float_t x_size, float_t y_size){
     cells_[index].x_size_ = x_size; cells_[index].y_size_ = y_size;
}
}; // namespace gp
}; // namespace coloquinte

#endif

