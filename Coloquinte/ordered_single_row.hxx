
#ifndef COLOQUINTE_OSRP
#define COLOQUINTE_OSRP

#include "common.hxx"

#include <queue>

namespace coloquinte{

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
    T get_displacement(legalizable_task<T> const newly_pushed, bool update){
        T target_abs_pos = newly_pushed.target_pos - current_width();
        T width = newly_pushed.width;
        T slope = - width;

        T cur_pos  = end;
        T cur_cost = 0;

        std::vector<OSRP_bound> passed_bounds;

        while( (slope < 0 or cur_pos > end - current_width())
           and not bounds.empty() and bounds.top().absolute_pos > target_abs_pos){
            slope += bounds.top().weight;
            T old_pos = cur_pos;
            cur_pos = bounds.top().absolute_pos;
            cur_cost += (old_pos - cur_pos) * (slope + width); // The additional cost for the other cells encountered

            // Remember which bounds we encountered in order to reset the object to its initial state
            if(not update)
                passed_bounds.push_back(bounds.top());
            bounds.pop();
        }

        T final_abs_pos = std::min(end - current_width() - width, // Always before the end
                                std::max(begin, slope >= 0 ? cur_pos : target_abs_pos) // Always after the beginning, but did we stop before reaching the target position? 
                                    );
        //assert(final_abs_pos >= begin);
        //assert(final_abs_pos <= end - current_width() - width);

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
    std::vector<result_t> get_placement() const{
        auto final_abs_pos = constraining_pos;
        std::partial_sum(final_abs_pos.rbegin(), final_abs_pos.rend(), final_abs_pos.rbegin(), [](T a, T b)->T{ return std::min(a,b); });

        std::vector<result_t> ret(cells.size());
        for(index_t i=0; i<cells.size(); ++i){
            ret[i] = result_t(cells[i], final_abs_pos[i] + prev_width[i]);
            //assert(final_abs_pos[i] >= begin);
            //assert(final_abs_pos[i] + prev_width[i+1] <= end);
        }
        return ret;
    }

};



} // End namespace coloquinte

#endif

