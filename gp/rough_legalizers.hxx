
#ifndef COLOQUINTE_GP_ROUGH_LEGALIZER
#define COLOQUINTE_GP_ROUGH_LEGALIZER

#include "common.hxx"

#include <vector>
#include <cassert>
#include <cmath>

/*
 * A simple class to perform rough legalization with extreme efficiency
 *
 *
 */

namespace coloquinte{
namespace gp{

class region_distribution{
    /*
     * Coordinates are mostly float but obstacles and areas are transformed to integers at construction for correctness
     */
    
    public:

    struct fixed_cell{
        box<int_t> box_;
        fixed_cell();
        fixed_cell(box<int_t> bx);
        /*
         * Later extension to handle whitespace with capacities different than area
         */
        // uint64_t capacity;
    };
    
    struct movable_cell{
        capacity_t demand_; // == area; No FP!!!
        float_t x_pos_, y_pos_;  // Target position, determining the cost to allocate it
        // int_t x_size, y_size; // May split cells
        index_t index_in_placement_;

        movable_cell(){}
        movable_cell(capacity_t demand, float_t x, float_t y, index_t ind) : demand_(demand), x_pos_(x), y_pos_(y), index_in_placement_(ind){}
    };

    private:

    struct region;
    struct cell_ref;
    
    struct cell_ref{
        capacity_t allocated_capacity_;
        float_t x_pos_, y_pos_;
        index_t index_in_list_;
        float_t marginal_cost_;

        cell_ref(){}
        cell_ref(capacity_t demand, float_t x, float_t y, index_t ind) : allocated_capacity_(demand), x_pos_(x), y_pos_(y), index_in_list_(ind){}

        bool operator<(cell_ref const o) const{ return marginal_cost_ < o.marginal_cost_; }
        friend region;
    };
    
    struct region{
        public:
        capacity_t capacity_, // ==area; No FP!!! 
            unused_capacity_;
        float_t x_pos_, y_pos_;
    
        box<int_t> surface_;
        std::vector<cell_ref> cell_references_;
        std::vector<fixed_cell> obstacles_;

        public:
        region(){} // Necessary if we want to resize vectors 
        region(box<int_t> bx, std::vector<fixed_cell> obstacles, std::vector<cell_ref> cells);
    
        void selfcheck() const;
        void x_bipartition(region & lft, region & rgt);
        void y_bipartition(region & up , region & dwn);

        float_t capacity() const;
        float_t allocated_capacity() const;
        float_t unused_capacity() const;
        float_t distance(point<float_t> const & P) const;
        float_t distance(cell_ref const & C) const;
        static void distribute_new_cells(region & a, region & b, std::vector<cell_ref> cells);
    };

    private:
    // Members
    //
    index_t x_cuts_cnt_, y_cuts_cnt_;
    
    box<int_t> placement_area_;
    std::vector<region> placement_regions_;
    std::vector<movable_cell> cell_list_;
    
    private:
    // Helper functions
    
    // Reduces the number of cuts in the current solution to at most region_cnt() - 1 without loss of solution quality
    void fractions_minimization();
    void redo_bipartition(region & Ra, region & Rb);
    std::vector<point<float_t> > get_exported_positions() const;

    void selfcheck() const;
    region & get_region(index_t x_coord, index_t y_coord);

    public:
    
    inline index_t x_regions_cnt() const;
    inline index_t y_regions_cnt() const;
    inline index_t regions_cnt()   const;
    
    // Dimensions of the uniform cells superimposed on the placement region
    inline float_t x_cell_dim() const;
    inline float_t y_cell_dim() const;
    
    /*
     * Two types of cost
     *    Region center estimation  : upper bound of legalization cost
     *    Export scaling estimation : lower bound of legalization cost
     */
    
    float_t cost() const;
    float_t exported_cost() const;
    
    /*
     * Further partitions
     */
    
    void x_bipartition();
    void y_bipartition();
    void quadpartition();
    
    /*
     * Optimization functions
     */
    
    // Improve bipartitions between closest non-empty neighbours
    void redo_bipartitions();
    
    // Tries to escape local minimas with long-distance moves to non-filled places
    void line_moves();
    
    /*
     * Create and export
     *
     * Create takes position and area infos from a global_placement object
     *
     * Export uses scaling of x/y coordinates
     *    Intelligent: tries to obtain minimal disruption
     *    But independent: illegal placement in a region is likely
     */
    
    //region_distribution(global_placement const & orig);
    //global_placement export_global_placement() const;
    
    region_distribution(box<int_t> placement_area, std::vector<movable_cell> all_cells);
};

inline region_distribution::fixed_cell::fixed_cell(){}
inline region_distribution::fixed_cell::fixed_cell(box<int_t> bx) : box_(bx){}

inline index_t region_distribution::x_regions_cnt() const { return 1 << x_cuts_cnt_; }
inline index_t region_distribution::y_regions_cnt() const { return 1 << y_cuts_cnt_; }
inline index_t region_distribution::regions_cnt()   const { index_t ret = x_regions_cnt() * y_regions_cnt(); assert(placement_regions_.size() == ret); return ret; }
inline region_distribution::region & region_distribution::get_region(index_t x_coord, index_t y_coord){
    assert(x_coord < x_regions_cnt() && y_coord < y_regions_cnt());
    return placement_regions_[y_coord * x_regions_cnt() + x_coord];
}

inline float_t region_distribution::region::capacity() const{ return capacity_; }
inline float_t region_distribution::region::unused_capacity() const{ return unused_capacity_; }
inline float_t region_distribution::region::allocated_capacity() const{
    capacity_t ret = 0;
    for(cell_ref const C : cell_references_){
       ret += C.allocated_capacity_; 
    }
    assert(unused_capacity() + ret == capacity());
    return ret;
}

inline float_t region_distribution::region::distance(point<float_t> const & P) const{
    float_t manhattan = std::abs(x_pos_ - P.x_) + std::abs(y_pos_ - P.y_);
    return manhattan * manhattan;
}
inline float_t region_distribution::region::distance(region_distribution::cell_ref const & C) const{
    return distance(point<float_t>(C.x_pos_, C.y_pos_));
}


} // Namespace gp
} // Namespace coloquinte

#endif

