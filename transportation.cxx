
#include "gp/transportation.hxx"


namespace coloquinte{
namespace gp{

void current_allocation::update_edge(index_t r1, index_t r2){
    while(not best_interregions_costs[r1][r2].empty()){
        movable_source cur = best_interregions_costs[r1][r2].top();
        // Test if the edge still exists
        if(sr_allocations[r1][cur.source] != 0){
            // Found the edge: stop
            break;
        }
        else{
            // This edge is in fact empty
            best_interregions_costs[r1][r2].pop();
        }
    }

    if(best_interregions_costs[r1][r2].empty()){
        // There is no edge: return
        return;
    }
    else{
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

void current_allocation::add_source_to_heaps(index_t r, index_t source){
    for(index_t i=0; i<region_cnt(); ++i){
        if(i == r) continue;
        best_interregions_costs[r][i].push(
            movable_source(source,
                sr_costs[i][source] - sr_costs[r][source]
            )
        );
    }
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

    assert(sr_allocations[reg][cur_source] >= 0); // The source to be pushed is indeed present in the region
    assert(r_capacities[reg] == 0); // The region is full, which explains why we need to push
    assert(flow <= arc_capacities[reg]); // The flow is not bigger than what can be sent

    // Only if the source was not already present here and we may have to push a source in the destination region
    if(not already_present and r_capacities[r_parents[reg]] == 0){ 
        // New source added to a region full region: rerun Dijkstra at the end
        add_source_to_heaps(r_parents[reg], cur_source);
        return true;
    }
    else if(arc_capacities[reg] == flow){
        // The source has been deleted from a region: rerun Dijkstra at the end
        return true;
    }
    else{
        // arc_capacities[reg] > flow
        arc_capacities[reg] -= flow;
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

capacity_t current_allocation::push_path(index_t pushed_reg, capacity_t demanded){
    // Get the final flow sent, which is smaller than the capacities on the path
    capacity_t flow = demanded;
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

    // If an edge changes cost or a region is full,
    // we need to update the costs, parents, sources and arc_capacities using a Dijkstra
    if(rerun_dijkstra)
        dijkstra_update();

    assert(flow > 0);
    return flow;
}

void current_allocation::add_source(capacity_t demand, std::vector<float_t> const & costs){
    index_t elt_ind = sr_allocations[0].size();

    for(index_t i=0; i<region_cnt(); ++i){
        assert(sr_costs[i].size() == elt_ind);
        assert(sr_allocations[i].size() == elt_ind);
        sr_costs[i].push_back(costs[i]);
        sr_allocations[i].push_back(0);
    }
    while(demand > 0){
        ++ dijkstra_cnt;
        index_t best_reg = null_ind;
        float_t best_cost = std::numeric_limits<float_t>::infinity();
        for(index_t reg=0; reg<region_cnt(); ++reg){
            // Find the region which gets the source
            if(r_costs[reg] + costs[reg] < best_cost){
                best_reg = reg;
                best_cost = r_costs[reg] + costs[reg];
            }
        }
        if(best_reg == null_ind){ throw std::runtime_error("No reachable region found\n"); }

        // Get the path's capacity and update the data structures
        capacity_t this_path_cap = push_path(best_reg, demand);

        // Substract the fulfilled demand
        demand -= this_path_cap;

        // Lazily store the change
        sr_allocations[best_reg][elt_ind] += this_path_cap;
    }

    // Set the source's demand
    for(index_t i=0; i<region_cnt(); ++i){
        if(r_capacities[i] == 0 and sr_allocations[i][elt_ind] > 0){
            add_source_to_heaps(i, elt_ind);
        }
    }
}

} // Namespace gp
} // Namespace coloquinte

