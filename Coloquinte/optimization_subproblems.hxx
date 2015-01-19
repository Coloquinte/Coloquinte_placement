
#ifndef COLOQUINTE_GP_OPTSUBPROBLEMS 
#define COLOQUINTE_GP_OPTSUBPROBLEMS 

#include "common.hxx"

#include <queue>
#include <vector>
#include <cassert>

namespace coloquinte{

typedef std::pair<float_t, capacity_t> t1D_elt;

std::vector<capacity_t>  optimize_1D(std::vector<t1D_elt> sources, std::vector<t1D_elt> sinks);

template<typename T>
struct legalizable_task{
    T width;
    T target_pos;
    index_t ind;
    legalizable_task(T w, T p, index_t i) : width(w), target_pos(p), ind(i){}
    bool operator<(legalizable_task<T> const o) const{ return target_pos < o.target_pos; }
};

// A class to obtain the optimal positions minimizing total weighted displacement along a row
// It is an ordered single row problem/fixed order single machine scheduling problem, solved by the clumping/specialized cascading descent algorithm
// The cost is linear in the distance to the target position, weighted by the width of the cells
template<typename T>
class OSRP_leg{
    struct OSRP_bound{
        T absolute_pos; // Will be the target absolute position of the cell
        T weight;       // Will be the width of the cell
    
        bool operator<(OSRP_bound const o) const{ return absolute_pos < o.absolute_pos; }
        OSRP_bound(T w, T abs_pos) : absolute_pos(abs_pos), weight(w) {}
    };

    T begin, end;

    std::vector<index_t> cells;            // The indexes in the circuit
    std::vector<T>   constraining_pos; // Where the cells have been pushed and constrain the positions of preceding cells
    std::vector<T>   prev_width;       // Cumulative width of the cells: calculates the absolute position of new cells

    std::priority_queue<OSRP_bound> bounds;

    // Get the cost of pushing a cell on the row
    T get_displacement(legalizable_task<T> const newly_pushed, bool update);

    public:
    T current_width() const{ return prev_width.back(); }
    T remaining_space() const{ return end - begin - current_width(); }
    T last_available_pos() const{ return constraining_pos.back() + current_width(); }

    T get_cost(legalizable_task<T> const task){ return get_displacement(task, false); }
    void push(legalizable_task<T> const task){ get_displacement(task, true); }

    // Initialize
    OSRP_leg(T b, T e) : begin(b), end(e), prev_width(1, 0) {}
    OSRP_leg(){}

    typedef std::pair<index_t, T> result_t;

    // Get the resulting placement
    std::vector<result_t> get_placement() const;
};


class full_single_row{
    struct bound{
        int_t abs_pos;
        float_t slope_diff;

        bool operator<(bound const o) const{ return abs_pos < o.abs_pos; }
        bound(int_t p, float_t s) : abs_pos(p), slope_diff(s) {}
    };

    float_t cur_slope;
    int_t lower, upper;

    std::vector<int_t> prev_width;
    std::vector<int_t> constraining_pos;
    std::priority_queue<bound> bounds;

    void update_positions();

    public:

    // Low-level functions to avoid internally building vectors
    void push_cell(int_t width, int_t lower_lim, int_t upper_lim);  // Give the characteristics for a cell
    void push_bound(int_t offset, float_t slope_diff);              // Push a bound for this cell
    void push_slope(float_t right_slope);                           // Push additional slope

    // Get the result
    std::vector<int_t> get_placement();

    full_single_row() : cur_slope(0.0), lower(std::numeric_limits<int_t>::min()), upper(std::numeric_limits<int_t>::max()), prev_width(1, 0) {}
};

template<typename T>
inline T OSRP_leg<T>::get_displacement(legalizable_task<T> const newly_pushed, bool update){
    T target_abs_pos = newly_pushed.target_pos - current_width();
    T width = newly_pushed.width;
    T slope = - width;

    T cur_pos  = end;
    T cur_cost = 0;

    std::vector<OSRP_bound> passed_bounds;

    while( not bounds.empty() and
        ((slope < 0 and bounds.top().absolute_pos > target_abs_pos) // Not reached equilibrium
        or bounds.top().absolute_pos > end - current_width() - width) // Still not a legal position
        ){
        T old_pos = cur_pos;
        cur_pos = bounds.top().absolute_pos;
        cur_cost += (old_pos - cur_pos) * (slope + width); // The additional cost for the other cells encountered
        slope += bounds.top().weight;

        // Remember which bounds we encountered in order to reset the object to its initial state
        if(not update)
            passed_bounds.push_back(bounds.top());
        bounds.pop();
    }

    T final_abs_pos = std::min(end - current_width() - width, // Always before the end and after the beginning
                            std::max(begin, slope >= 0 ? cur_pos : target_abs_pos) // but did we stop before reaching the target position? 
                                );

    cur_cost += (cur_pos - final_abs_pos) * (slope + width); // The additional cost for the other cells encountered

    if(std::numeric_limits<T>::is_integer){
        assert(final_abs_pos >= begin);
        assert(final_abs_pos <= end - current_width() - width);
    }

    if(update){
        prev_width.push_back(width + current_width());
        cells.push_back(newly_pushed.ind);
        constraining_pos.push_back(final_abs_pos);
        if(slope > 0){ // Remaining capacity of an encountered bound
            bounds.push(OSRP_bound(slope, cur_pos));
        }
        // The new bound, minus what it absorbs of the remaining slope
        if(target_abs_pos > begin){
            bounds.push(OSRP_bound(2*width + std::min(slope, static_cast<T>(0) ), target_abs_pos));
        }
    }
    else{
        for(OSRP_bound b : passed_bounds){
            bounds.push(b);
        }
    }

    return cur_cost + width * std::abs(final_abs_pos - target_abs_pos); // Add the cost of the new cell
}


inline void full_single_row::update_positions(){
    int_t cur_pos = upper;
    // If we didn't push the position of the row
    if(constraining_pos.size() + 1 < prev_width.size()){
        while(not bounds.empty() and (cur_slope > 0.0 or bounds.top().abs_pos > upper)){
            cur_slope -= bounds.top().slope_diff;
            cur_pos = bounds.top().abs_pos;
            bounds.pop();
        }
        int_t final_abs_pos = std::max(std::min(cur_pos, upper), lower);
        constraining_pos.push_back(final_abs_pos);
        if(cur_slope < 0.0){
            bounds.push(bound(final_abs_pos, -cur_slope));
        }
    }
}

inline void full_single_row::push_cell(int_t width, int_t lower_lim, int_t upper_lim){
    update_positions();

    lower = std::max(lower, lower_lim - prev_width.back());
    prev_width.push_back(width + prev_width.back());
    upper = upper_lim - prev_width.back();
    cur_slope = 0.0;
}

inline void full_single_row::push_slope(float_t right_slope){
    cur_slope += right_slope;
}

inline void full_single_row::push_bound(int_t position, float_t slope_diff){
    assert(constraining_pos.size() + 1 < prev_width.size());
    bounds.push(bound(position - prev_width[prev_width.size()-2], slope_diff));
}

inline std::vector<int_t> full_single_row::get_placement(){
    update_positions();

    auto vals = std::vector<int_t>(constraining_pos.size());
    std::partial_sum(constraining_pos.rbegin(), constraining_pos.rend(), vals.rbegin(), [](int_t a, int_t b)->int_t{ return std::min(a,b); });
    for(index_t i=0; i<vals.size(); ++i){
        vals[i] += prev_width[i];
    }
    return vals;
}

template<typename T>
inline std::vector<std::pair<index_t, T> > OSRP_leg<T>::get_placement() const{
    auto final_abs_pos = constraining_pos;
    std::partial_sum(final_abs_pos.rbegin(), final_abs_pos.rend(), final_abs_pos.rbegin(), [](T a, T b)->T{ return std::min(a,b); });

    std::vector<result_t> ret(cells.size());
    for(index_t i=0; i<cells.size(); ++i){
        ret[i] = result_t(cells[i], final_abs_pos[i] + prev_width[i]);

        if(std::numeric_limits<T>::is_integer){
            assert(final_abs_pos[i] >= begin);
            assert(final_abs_pos[i] + prev_width[i+1] <= end);
        }
    }
    return ret;
}

class current_allocation{
    static const index_t null_ind = std::numeric_limits<index_t>::max();

    // Internal data structures

    // Priority queue element to determine the source to be used between regions
    struct movable_source{
        index_t source;
        float_t cost;
        bool operator<(movable_source const o) const{
               return cost > o.cost // Sorted by cost
            || (cost == o.cost && source < o.source); // And by index to limit the number of fractional elements between two regions
        }
        movable_source(index_t s, float_t c) : source(s), cost(c) {}
    };

    // Member data


    // The current state
    std::vector<std::vector<capacity_t>  > sr_allocations; // For each region, for each source, the capacity allocated by the region
    std::vector<std::vector<float_t> >     sr_costs;       // The costs from a region to a source
    std::vector<capacity_t>                r_capacities;   // The remaining capacities of the regions

    // Shortest path data
    std::vector<float_t>                   r_costs;        // The costs of allocating to a region
    std::vector<index_t>                   r_parents;      // The parents of the regions i.e. the regions where we push sources first (or null_ind)
    std::vector<index_t>                   r_sources;      // The sources involved in these edges
    std::vector<capacity_t>                arc_capacities; // The capacities of the edges to the parents, or of the region if no parent

    // Best edges data
    std::vector<std::vector<std::priority_queue<movable_source> > > best_interregions_costs; // What is the best source to move to go from region k1 to region k2?
    index_t dijkstra_cnt;


    // Helper functions

    // Number of regions
    index_t region_cnt() const{
        assert(sr_costs.size() == sr_allocations.size());
        return sr_costs.size();
    }

    // Update the edge between two regions
    void update_edge(index_t r1, index_t r2);
    // Add a source to all heaps of a region; returns if we need to update a path
    bool add_source_to_heaps(index_t r, index_t source);
    // Initialize the heaps of a region
    void create_heaps(index_t reg);

    // Run the shortest path algorithm to update the cost of each region
    void dijkstra_update();

    // Update the edge and returns if we need to rerun Dijkstra
    bool push_edge(index_t reg, capacity_t flow);
    // Updates a full path when pushing an element; returns if we need to rerun Dijkstra
    bool push_path(index_t pushed_reg, capacity_t demanded, capacity_t & flow);

    public:
    // Add a new source to the transportation problem; should be done in decreasing order of demand to keep low complexity
    void add_source(capacity_t demand, std::vector<float_t> const & costs);

    current_allocation(std::vector<capacity_t> caps)
        :
        sr_allocations(caps.size()),
        sr_costs(caps.size()),
        r_capacities(caps),
        r_costs(caps.size(), 0.0),
        r_parents(caps.size(), null_ind),
        r_sources(caps.size(), null_ind),
        arc_capacities(caps),
        best_interregions_costs(caps.size(), std::vector<std::priority_queue<movable_source> >(caps.size())),
        dijkstra_cnt(0)
        {
            assert(caps.size() > 0);
            dijkstra_update();
        }

    std::vector<std::vector<capacity_t> > get_allocations() const{ return sr_allocations; }
    index_t get_iterations_cnt() const { return dijkstra_cnt; }
};

#endif

}

