
#include "Coloquinte/optimization_subproblems.hxx"


namespace coloquinte{

std::vector<capacity_t>  transport_1D(std::vector<t1D_elt> sources, std::vector<t1D_elt> sinks){
    /* Description of the algorithm:
     *
     *    For each cell, put it in its optimal region or the last region where a cell is if there is no space
     *    Push the changes in the derivative of the cost function to a priority queue; those changes occur
     *          when evicting the preceding cell from its current region
     *          when moving to a non-full region
     *    While the new cell would occupy region which was still free, get the new slope (derivative)
     *    and push all preceding cell until this region is freed or the slope is 0
     */

    struct bound{
        capacity_t pos;
        float_t slope_diff;
        bool operator<(bound const o) const{ return pos < o.pos; }
    };

    std::priority_queue<bound> bounds;
    std::vector<capacity_t> constraining_pos;
    std::vector<capacity_t> prev_cap(1, 0), prev_dem(1, 0);
    for(auto const s : sinks){
        prev_cap.push_back(s.second + prev_cap.back());
    }
    for(auto const s : sources){
        prev_dem.push_back(s.second + prev_dem.back());
    }
    assert(prev_cap.back() >= prev_dem.back());

    const capacity_t min_abs_pos = 0, max_abs_pos = prev_cap.back() - prev_dem.back();
    assert(min_abs_pos <= max_abs_pos);

    auto push_bound = [&](capacity_t p, float_t s){
        assert(s >= -0.0);
        if(p > min_abs_pos){
            bound B;
            B.pos = p;
            B.slope_diff = s;
            bounds.push(B);
        }
    }; 

    // Distance to the right - distance to the left
    auto get_slope = [&](index_t src, index_t boundary){
        return std::abs(sources[src].first - sinks[boundary+1].first) - std::abs(sources[src].first - sinks[boundary].first);
    };

    capacity_t cur_abs_pos = min_abs_pos;
    index_t opt_r=0, next_r=0, first_free_r=0;

    for(index_t i=0; i<sources.size(); ++i){
        // Update the optimal region
        while(opt_r+1 < sinks.size() and 0.5 * (sinks[opt_r].first + sinks[opt_r+1].first) < sources[i].first){
            ++opt_r;
        }
        // Update the next region
        index_t prev_next_r = next_r;
        while(next_r < sinks.size() and sinks[next_r].first < sources[i].first){
            ++next_r;
        }

        if(i>0){
            //for(index_t j=std::max(prev_next_r,1u)-1; j<std::min(first_free_r, next_r+1); ++j){
            for(index_t j=std::max(prev_next_r,1u)-1; j<std::min(first_free_r, opt_r); ++j){
                assert(get_slope(i,j) <= get_slope(i-1,j));
                push_bound(prev_cap[j+1] - prev_dem[i], get_slope(i-1, j) - get_slope(i,j));
            }
        }
        // Add the bounds due to crossing the boundaries alone
        for(index_t j=first_free_r; j<opt_r; ++j){
            assert(get_slope(i,j) <= 0.0);
            push_bound(prev_cap[j+1] - prev_dem[i], -get_slope(i, j));
        }

        first_free_r = std::max(first_free_r, opt_r);
        capacity_t this_abs_pos = std::max(cur_abs_pos, prev_cap[first_free_r] - prev_dem[i]); // Just after the previous cell or at the beginning of the destination region

        while(first_free_r+1 < sinks.size() and this_abs_pos > std::max(prev_cap[first_free_r+1] - prev_dem[i+1], min_abs_pos)){ // Absolute position that wouldn't make the cell fit in the region, and we are not in the last region yet
            capacity_t end_pos = std::max(prev_cap[first_free_r+1] - prev_dem[i+1], min_abs_pos);

            float_t add_slope = get_slope(i, first_free_r);
            float_t slope = add_slope;

            while(not bounds.empty() and slope >= 0.0 and bounds.top().pos > end_pos){
                this_abs_pos = bounds.top().pos;
                slope -= bounds.top().slope_diff;
                bounds.pop();
            }
            if(slope >= 0.0){ // We still push: the cell completely escapes the region
                this_abs_pos = end_pos;
                push_bound(end_pos, add_slope-slope);
            }
            else{ // Ok, absorbed the whole slope: push what remains and we still occupy the next region
                push_bound(this_abs_pos, -slope);
                ++first_free_r;
            }
        }
        cur_abs_pos = this_abs_pos;
        constraining_pos.push_back(this_abs_pos);
    }

    assert(constraining_pos.size() == sources.size());

    if(not constraining_pos.empty()){
        // Calculate the final constraining_pos
        constraining_pos.back() = std::min(max_abs_pos, constraining_pos.back());
    }

    std::partial_sum(constraining_pos.rbegin(), constraining_pos.rend(), constraining_pos.rbegin(), [](capacity_t a, capacity_t b)->capacity_t{ return std::min(a, b); });

    for(index_t i=0; i<constraining_pos.size(); ++i){
        constraining_pos[i] += prev_dem[i];
    }

    return constraining_pos;
}

namespace{ // Anonymous namespace to hide the transportation structures

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
    std::vector<capacity_t>                s_demands;      // The demands of the sources
    std::vector<capacity_t>                r_capacities;   // The remaining capacities of the regions

    // Shortest path data
    std::vector<float_t>                   r_costs;        // The costs of allocating to a region
    std::vector<index_t>                   r_parents;      // The parents of the regions i.e. the regions where we push sources first (or null_ind)
    std::vector<index_t>                   r_sources;      // The source involved in these edges
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
    void add_source(index_t elt_ind);

    current_allocation(std::vector<capacity_t> caps, std::vector<capacity_t> demands, std::vector<std::vector<float_t> > costs)
        :
        sr_allocations(caps.size()),
        sr_costs(costs),
        s_demands(demands),
        r_capacities(caps),
        r_costs(caps.size(), 0.0),
        r_parents(caps.size(), null_ind),
        r_sources(caps.size(), null_ind),
        arc_capacities(caps),
        best_interregions_costs(caps.size(), std::vector<std::priority_queue<movable_source> >(caps.size())),
        dijkstra_cnt(0)
        {
            assert(caps.size() > 0);
            assert(costs.size() == caps.size());
            dijkstra_update();
        }

    std::vector<std::vector<capacity_t> > get_allocations() const{ return sr_allocations; }
    index_t get_iterations_cnt() const { return dijkstra_cnt; }
};

void current_allocation::update_edge(index_t r1, index_t r2){
    while(not best_interregions_costs[r1][r2].empty() and sr_allocations[r1][best_interregions_costs[r1][r2].top().source] == 0){
        best_interregions_costs[r1][r2].pop();
    }

    if(not best_interregions_costs[r1][r2].empty()){
        // There is an edge
        movable_source cur = best_interregions_costs[r1][r2].top();
        float_t new_cost = r_costs[r2] + cur.cost;
        if(new_cost < r_costs[r1]){
            r_costs[r1] = cur.cost;
            r_sources[r1] = cur.source;
            r_parents[r1] = r2;
            arc_capacities[r1] = sr_allocations[r1][cur.source];
        }
    }
}

bool current_allocation::add_source_to_heaps(index_t r, index_t source){
    bool need_rerun = false;
    for(index_t i=0; i<region_cnt(); ++i){
        if(i == r) continue;
        best_interregions_costs[r][i].push(
            movable_source(source,
                sr_costs[i][source] - sr_costs[r][source]
            )
        );
        while(sr_allocations[r][best_interregions_costs[r][i].top().source] == 0){
            best_interregions_costs[r][i].pop();
        }
        need_rerun = (best_interregions_costs[r][i].top().source == source) or need_rerun;
    }
    return need_rerun;
}

void current_allocation::create_heaps(index_t reg){
    // Get all relevant elements
    std::vector<std::vector<movable_source> > interregion_costs(region_cnt());
    for(index_t i=0; i<sr_allocations[reg].size(); ++i){
        if(sr_allocations[reg][i] > 0){
            for(index_t oreg=0; oreg<region_cnt(); ++oreg){
                if(oreg == reg) continue;
                interregion_costs[oreg].push_back(
                    movable_source(
                        i,
                        sr_costs[oreg][i] - sr_costs[reg][i]
                    )
                );
            }
        }
    }
    // Create the heaps
    for(index_t oreg=0; oreg<region_cnt(); ++oreg){
        best_interregions_costs[reg][oreg] = std::priority_queue<movable_source>(interregion_costs[oreg].begin(), interregion_costs[oreg].end());
    }
}

// Returns if the path has been modified so that we would need to rerun Dijkstra
bool current_allocation::push_edge(index_t reg, capacity_t flow){
    index_t cur_source = r_sources[reg];

    // Does this edge allocates a new source in the destination region? If yes, update the corresponding heaps
    bool already_present = sr_allocations[r_parents[reg]][cur_source] > 0;

    // Deallocating from the first region is handled by the get_edge function: just substract the flow
    sr_allocations[          reg ][cur_source] -= flow;
    sr_allocations[r_parents[reg]][cur_source] += flow;

    assert(sr_allocations[reg][cur_source] >= 0); // The source to be pushed was indeed present in the region
    assert(r_capacities[reg] == 0); // The region is full, which explains why we need to push
    assert(flow <= arc_capacities[reg]); // The flow is not bigger than what can be sent

    arc_capacities[reg] = sr_allocations[reg][cur_source]; // Just update the capacity if it turns out that we don't need to run Dijkstra
    
    if(arc_capacities[reg] == 0){
        // The source may have been deleted from a region: rerun Dijkstra at the end
        return true;
    }
    else if(not already_present and r_capacities[r_parents[reg]] == 0){
        // A new source is allocated to a full region: rerun Dijkstra at the end if it changed the heap's top
        return add_source_to_heaps(r_parents[reg], cur_source);
    }
    else{
        // The edge is still present with the same cost and non-zero updated capacity
        // The path still exists: no need to rerun Dijkstra yet
        return false;
    }
}

void current_allocation::dijkstra_update(){
    // Simple case of the regions with remaining capacity
    std::vector<int> visited(region_cnt(), 0);
    index_t visited_cnt = 0;
    for(index_t i=0; i<region_cnt(); ++i){
        r_sources[i] = null_ind;
        r_parents[i] = null_ind;
        if(r_capacities[i] > 0){
            r_costs[i] = 0.0;
            arc_capacities[i] = r_capacities[i];

            visited[i] = 1;
            ++visited_cnt;
        }
        else{
            r_costs[i] = std::numeric_limits<float_t>::infinity();
            arc_capacities[i] = 0;
        }
    }
    // if(visited_cnt <= 0) throw std::runtime_error("Capacity problem: no region has been marked as reachable\n");
    if(visited_cnt == region_cnt()){ return; }
    // Get the costs for every non-visited region
    for(index_t i=0; i<region_cnt(); ++i) if(visited[i] == 0){ // For every region that is not visited yet
        for(index_t j=0; j<region_cnt(); ++j) if(visited[j] == 1){ // For every already visited region
            // Get the best interregion cost
            update_edge(i,j);
        }
    }
    while(visited_cnt < region_cnt()){
        // Find the region with the lowest cost to visit; mark it visited
        index_t best_reg = null_ind;
        float_t best_cost = std::numeric_limits<float_t>::infinity();
        for(index_t i=0; i<region_cnt(); ++i) if(visited[i] == 0){ // For every region that is not visited yet
            if(r_costs[i] < best_cost){
                best_cost = r_costs[i];
                best_reg  = i;
            }
        }
        if(best_reg == null_ind) break; // Some regions are unreachable, typically because they have zero capacity at the beginning
        visited[best_reg] = 1;
        ++visited_cnt;
        // Update the cost for every unvisited region
        for(index_t i=0; i<region_cnt(); ++i) if(visited[i] == 0){ // For every region that is not visited yet
            update_edge(i, best_reg);
        }
    }
}

bool current_allocation::push_path(index_t pushed_reg, capacity_t demanded, capacity_t & flow){
    // Get the final flow sent, which is smaller than the capacities on the path
    flow = demanded;
    for(index_t reg = pushed_reg; reg != null_ind; reg = r_parents[reg]){
        flow = std::min(flow, arc_capacities[reg]);
    }

    bool rerun_dijkstra = false;
    // Update the path between the regions
    index_t reg = pushed_reg;
    for(; r_parents[reg] != null_ind; reg = r_parents[reg]){
        assert(r_capacities[reg] == 0);
        rerun_dijkstra = push_edge(reg, flow) or rerun_dijkstra;
    }

    assert(r_capacities[reg] > 0);
    assert(arc_capacities[reg] == r_capacities[reg]);
    assert(r_capacities[reg] >= flow);

    // Update the capacities at the end
    r_capacities[reg] -= flow;
    arc_capacities[reg] -= flow;

    // The last region on the path is the one that satisfies the demand
    if(r_capacities[reg] == 0){ // If we just consumed the available capacity, it becomes useful to move sources off this region: build the heap
        create_heaps(reg);
        rerun_dijkstra = true;
    }

    assert(flow > 0);

    // If an edge changes cost or a region is full,
    // we need to update the costs, parents, sources and arc_capacities using a Dijkstra
    // but later
    return rerun_dijkstra;
}

void current_allocation::add_source(index_t elt_ind){ //capacity_t demand, std::vector<float_t> const & costs){
    for(index_t i=0; i<region_cnt(); ++i){
        sr_allocations[i].push_back(0);
    }

    bool need_rerun = false;
    capacity_t demand = s_demands[elt_ind];

    while(demand > 0){
        // In case we modified the structures earlier
        if(need_rerun){
            dijkstra_update();
            need_rerun = false;
        }

        ++ dijkstra_cnt;
        index_t best_reg = null_ind;
        float_t best_cost = std::numeric_limits<float_t>::infinity();
        for(index_t reg=0; reg<region_cnt(); ++reg){
            // Find the region which gets the source
            if(r_costs[reg] + sr_costs[reg][elt_ind] < best_cost){
                best_reg = reg;
                best_cost = r_costs[reg] + sr_costs[reg][elt_ind];
            }
        }
        if(best_reg == null_ind){ throw std::runtime_error("No reachable region found\n"); }

        capacity_t flow = 0;
        // Tells us whether we need to update the data structures
        need_rerun = push_path(best_reg, demand, flow);
        demand -= flow;

        // Lazily store the change
        sr_allocations[best_reg][elt_ind] += flow;
    }

    // Set the source's demand
    for(index_t i=0; i<region_cnt(); ++i){
        if(r_capacities[i] == 0 and sr_allocations[i][elt_ind] > 0){
            need_rerun = add_source_to_heaps(i, elt_ind) or need_rerun;
        }
    }
    // We leave a clean set with correct paths for the next iteration
    if(need_rerun)
        dijkstra_update();
}

} // End anonymous namespace

std::vector<std::vector<capacity_t> > transport_generic(std::vector<capacity_t> const & capacities, std::vector<capacity_t> const & demands, std::vector<std::vector<float_t> > const & costs){
    current_allocation transporter(capacities, demands, costs);

    for(index_t i=0; i<demands.size(); ++i){
        transporter.add_source(i);
    }

    return transporter.get_allocations();
}

} // Namespace coloquinte

