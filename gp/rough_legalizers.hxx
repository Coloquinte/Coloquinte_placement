
#ifndef COLOQUINTE_GP_ROUGH_LEGALIZER
#define COLOQUINTE_GP_ROUGH_LEGALIZER

#include "placement.hxx"

#include <vector>

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
    };
    
    struct cell_ref{
        capacity_t allocated_capacity_;
        float_t x_pos_, y_pos_;
        index_t index_in_list_;
        float_t marginal_cost_;
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
        region(box<int_t> bx, std::vector<fixed_cell> obstacles, std::vector<cell_ref> cells);
    
        void selfcheck() const;
        void x_bipartition(region & lft, region & rgt);
        void y_bipartition(region & up , region & dwn);

        private:
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
    void redo_bipartition(index_t indexa, index_t indexb);
    void selfcheck() const;
    
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
    float_t export_cost() const;
    
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
    
    // Move cells to the closest non-filled neighbours
    void neighbour_moves();
    
    // Tries to escape local minimas with long-distance moves
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
    
    region_distribution(global_placement const & orig);
    global_placement export_global_placement() const;
    
};

inline region_distribution::fixed_cell::fixed_cell(){}
inline region_distribution::fixed_cell::fixed_cell(box<int_t> bx) : box_(bx){}

inline index_t region_distribution::x_regions_cnt() const { return 1 << x_cuts_cnt_; }
inline index_t region_distribution::y_regions_cnt() const { return 1 << y_cuts_cnt_; }
inline index_t region_distribution::regions_cnt()   const { return x_regions_cnt() * y_regions_cnt(); }

} // Namespace gp
} // Namespace coloquinte

#endif

