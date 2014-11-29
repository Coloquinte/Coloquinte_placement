
#include "gp/rough_legalizers.hxx"
#include <algorithm>
#include <cmath>
#include <cassert>
#include <queue>

namespace coloquinte{
namespace gp{

const index_t null_ind = std::numeric_limits<index_t>::max();

class current_allocation{
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
    void update_edge(index_t r1, index_t r2){
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

    void add_source_to_heaps(index_t r, index_t source){
        for(index_t i=0; i<region_cnt(); ++i){
            if(i == r) continue;
            best_interregions_costs[r][i].push(
                movable_source(source,
                    sr_costs[i][source] - sr_costs[r][source]
                )
            );
        }
    }

    void create_heaps(index_t reg){
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
    bool push_edge(index_t reg, capacity_t flow){
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

    void dijkstra_update(){
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

    capacity_t push_path(index_t pushed_reg, capacity_t demanded){
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


    public:
    // Better: keep the cost of each region and its parent; update only what is changed: equivelent to inverting the Dijkstra algorithm, from regions to source
    // Nothing needs to be done until a region has no capacity left; when it happens, we need to create the heaps for this region
    // When no edge is destroyed, no changes are made; if an edge is destroyed (i.e. a source completely changes ownership), we need to update the costs of the zero-capacity regions
    void add_source(capacity_t demand, std::vector<float_t> const & costs){
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


void region_distribution::region::selfcheck() const{
    capacity_t total_allocated = 0;
    for(cell_ref const c : cell_references_){
        total_allocated += c.allocated_capacity_;
        if(c.allocated_capacity_ <= 0){ abort(); }
    }
    if(total_allocated > capacity_){ abort(); }
    for(index_t i=0; i+1<obstacles_.size(); ++i){
        for(index_t j=i+1; j<obstacles_.size(); ++j){
            assert(not obstacles_[i].box_.intersects(obstacles_[j].box_));
        }
    }
}

void region_distribution::region::uniquify_references(){
    std::sort(cell_references_.begin(), cell_references_.end(), [](cell_ref a, cell_ref b){ return a.index_in_list_ < b.index_in_list_; });

    std::vector<cell_ref> new_refs;
    if(cell_references_.size() >= 1){
        cell_ref prev_cell = cell_references_[0];
        for(auto it = cell_references_.begin()+1; it != cell_references_.end(); ++it){
            if(it->index_in_list_ == prev_cell.index_in_list_){
                prev_cell.allocated_capacity_ += it->allocated_capacity_;
            }
            else{
                new_refs.push_back(prev_cell);
                prev_cell = *it;
            }
        }
        new_refs.push_back(prev_cell);
    }
    std::swap(cell_references_, new_refs);
}

void region_distribution::selfcheck() const{
    for(region const & R : placement_regions_){
        R.selfcheck();
    }
    std::vector<capacity_t> capacities(cell_list_.size(), 0);
    for(region const & R : placement_regions_){
        for(cell_ref const C : R.cell_references_){
            capacities[C.index_in_list_] += C.allocated_capacity_;
            if(C.allocated_capacity_ <= 0){ abort(); }
        }
    }
    for(index_t i=0; i < cell_list_.size(); ++i){
        if(capacities[i] != cell_list_[i].demand_){ abort(); }
    }
}

region_distribution::region::region(box<int_t> bx, std::vector<fixed_cell> obstacles, std::vector<cell_ref> cells) : surface_(bx), cell_references_(cells){
    pos_ = static_cast<float_t>(0.5) * static_cast<point<float_t> >(
        point<int_t>(surface_.x_max_ + surface_.x_min_, surface_.y_max_ + surface_.y_min_)
    );

    capacity_ = static_cast<capacity_t>(surface_.x_max_ - surface_.x_min_) * static_cast<capacity_t>(surface_.y_max_ - surface_.y_min_);
    for(auto const & O : obstacles){
        box<int_t> const B = O.box_;
        if(surface_.intersects(B)){
            box<int_t> common_surface = surface_.intersection(B);
            obstacles_.push_back(common_surface);
            capacity_ -= static_cast<capacity_t>(common_surface.x_max_ - common_surface.x_min_) * static_cast<capacity_t>(common_surface.y_max_ - common_surface.y_min_);
        }
    }
}

void region_distribution::region::x_bipartition(region & lft, region & rgt){
    int_t x_min = surface_.x_min_,
          y_min = surface_.y_min_,
          x_max = surface_.x_max_,
          y_max = surface_.y_max_;
    int_t middle = (x_min + x_max) / 2; // x_max - x_min >= 0
    box<int_t> bx_lft(x_min, middle, y_min, y_max), bx_rgt(middle, x_max, y_min, y_max);
    
    lft = region(bx_lft, obstacles_, std::vector<cell_ref>());
    rgt = region(bx_rgt, obstacles_, std::vector<cell_ref>());

    distribute_new_cells(lft, rgt, cell_references_);

    assert(lft.allocated_capacity() + rgt.allocated_capacity() == allocated_capacity());
    assert(lft.capacity() + rgt.capacity() == capacity());
}

void region_distribution::region::y_bipartition(region & up, region & dwn){
    int_t x_min = surface_.x_min_,
          y_min = surface_.y_min_,
          x_max = surface_.x_max_,
          y_max = surface_.y_max_;
    int_t middle = (y_min + y_max) / 2; // x_max - x_min >= 0
    box<int_t> bx_up(x_min, x_max, middle, y_max), bx_dwn(x_min, x_max, y_min, middle);
    
    dwn = region(bx_dwn, obstacles_, std::vector<cell_ref>());
    up  = region(bx_up , obstacles_, std::vector<cell_ref>());

    distribute_new_cells(up, dwn, cell_references_);

    assert(up.allocated_capacity() + dwn.allocated_capacity() == allocated_capacity());
    assert(up.capacity() + dwn.capacity() == capacity());
}


// The big awful function that handles optimal cell distribution between two regions; not meant to be called externally
void region_distribution::region::distribute_new_cells(region & region_a, region & region_b, std::vector<cell_ref> cells){

    for(cell_ref & c : cells){
        c.marginal_cost_ = region_a.distance(c) - region_b.distance(c);
    }

    // Cells trending toward a first
    std::sort(cells.begin(), cells.end());//, [](cell_ref const c1, cell_ref const c2) -> bool{ return c1.marginal_cost_ < c2.marginal_cost_; });

    index_t preference_limit=0,         // First cell that would rather go to b (or cells.size())
         a_capacity_limit=0,            // After the last cell that region_a can take entirely (or 0)
         b_capacity_limit=cells.size(); // Last cell (but first in the vector) that region_b can take entirely (or cells.size())

    capacity_t remaining_capacity_a = region_a.capacity_, remaining_capacity_b = region_b.capacity_;
    for(;preference_limit < cells.size() && cells[preference_limit].marginal_cost_ <= 0.0; ++preference_limit);

    { // Block
    capacity_t remaining_cap_a = region_a.capacity_;
    index_t i=0;
    while(i<cells.size()){
        remaining_cap_a -= cells[i].allocated_capacity_;
        if(remaining_cap_a >= 0){
            a_capacity_limit = i+1;
            remaining_capacity_a = remaining_cap_a;
        }
        ++i;
    }
    } // Block

    { // Block
    capacity_t remaining_cap_b = region_b.capacity_;
    index_t i=cells.size();
    while(i>0){
        --i;
        remaining_cap_b -= cells[i].allocated_capacity_;
        if(remaining_cap_b >= 0){
            b_capacity_limit = i;
            remaining_capacity_b = remaining_cap_b;
        }
    }
    } // Block

    std::vector<cell_ref> & cells_a_side = region_a.cell_references_;
    std::vector<cell_ref> & cells_b_side = region_b.cell_references_;
    cells_a_side.clear();
    cells_b_side.clear();
    if(preference_limit >= b_capacity_limit and preference_limit <= a_capacity_limit){
        cells_a_side.insert(cells_a_side.end(), cells.begin(), cells.begin() + preference_limit);
        cells_b_side.insert(cells_b_side.end(), cells.begin() + preference_limit, cells.end());
    }
    else{
        index_t cut_position;
        capacity_t allocated_to_a_part;

        // Where do we cut?
        if(preference_limit < b_capacity_limit){ // Pack on b
            cut_position = b_capacity_limit-1; // Exists since preference_limit >= 0
            allocated_to_a_part = cells[cut_position].allocated_capacity_ - remaining_capacity_b;
        }
        else{ // Pack on a
        // if(preference_limit > a_capacity_limit)
            cut_position = a_capacity_limit; // Necessarily a correct position since preference_limits <= cells.size()
            allocated_to_a_part = remaining_capacity_a;
        }

        cells_a_side.reserve( cut_position             +1);
        cells_b_side.reserve(cells.size()-cut_position +1);

        cells_a_side.insert(cells_a_side.end(), cells.begin(), cells.begin() + cut_position);

        cell_ref cell_cut_a = cells[cut_position], cell_cut_b = cells[cut_position];
        cell_cut_a.allocated_capacity_ = allocated_to_a_part;
        cell_cut_b.allocated_capacity_ -= allocated_to_a_part;
        if(cell_cut_a.allocated_capacity_ > 0){ cells_a_side.push_back(cell_cut_a); }
        if(cell_cut_b.allocated_capacity_ > 0){ cells_b_side.push_back(cell_cut_b); }

        cells_b_side.insert(cells_b_side.end(), cells.begin() + cut_position+1, cells.end());
    }
}

void region_distribution::x_bipartition(){
    std::vector<region> old_placement_regions(2*regions_cnt());
    placement_regions_.swap(old_placement_regions);

    index_t old_x_regions_cnt = x_regions_cnt();
    index_t old_y_regions_cnt = y_regions_cnt();
    x_cuts_cnt_++;

    for(index_t x=0; x < old_x_regions_cnt; ++x){
        for(index_t y=0; y < old_y_regions_cnt; ++y){
            index_t i = y * old_x_regions_cnt + x;
            old_placement_regions[i].x_bipartition(get_region(2*x, y), get_region(2*x+1, y));
        }
    }
}

void region_distribution::y_bipartition(){
    std::vector<region> old_placement_regions(2*regions_cnt());
    placement_regions_.swap(old_placement_regions);

    index_t old_x_regions_cnt = x_regions_cnt();
    index_t old_y_regions_cnt = y_regions_cnt();
    y_cuts_cnt_++;

    for(index_t x=0; x < old_x_regions_cnt; ++x){
        for(index_t y=0; y < old_y_regions_cnt; ++y){
            index_t i = y * old_x_regions_cnt + x;
            old_placement_regions[i].y_bipartition(get_region(x, 2*y), get_region(x, 2*y+1));
        }
    }
}


void region_distribution::quadpartition(){
    // Naive implementation; maybe advanced ones could use an exact solution
    x_bipartition();
    y_bipartition();
}

void region_distribution::redo_bipartition(region_distribution::region & Ra, region_distribution::region & Rb){
    std::vector<cell_ref> cells;
    cells.reserve(Ra.cell_references_.size()+Rb.cell_references_.size());
    cells.insert(cells.end(), Ra.cell_references_.begin(), Ra.cell_references_.end());
    cells.insert(cells.end(), Rb.cell_references_.begin(), Rb.cell_references_.end());

    region::distribute_new_cells(Ra, Rb, cells);
}

void region_distribution::redo_bipartitions(){
    // This function performs optimization between neighbouring regions in various directions
    // The most important feature is diagonal optimization, since it is not done during partitioning
    // In order to optimize past obstacles even if only local optimization is considered, regions with no capacity are ignored

    // Perform optimization on the diagonals going South-East
    for(index_t diag = 0; diag < x_regions_cnt() + y_regions_cnt() - 1; ++diag){
        index_t x_begin, y_begin, nbr_elts;
        if(diag < y_regions_cnt()){ // Begin on the left side
            x_begin = 0;
            y_begin = diag;
            nbr_elts = std::min(x_regions_cnt(), y_regions_cnt() - diag);
        }
        else{ // Begin on the upper side
            x_begin = diag - y_regions_cnt() + 1;
            y_begin = 0;
            nbr_elts = std::min(y_regions_cnt(), x_regions_cnt() + y_regions_cnt() - diag -1);
        }
        index_t offs = 0;
        for(index_t nxt_offs = 1; nxt_offs < nbr_elts; ++nxt_offs){
            if(get_region(x_begin+nxt_offs, y_begin+nxt_offs).capacity() > 0){
                redo_bipartition(get_region(x_begin+offs, y_begin+offs), get_region(x_begin+nxt_offs, y_begin+nxt_offs));
                offs = nxt_offs;
            }
        }
    }

    // Perform optimization on the diagonals going South-West
    for(index_t diag = 0; diag < x_regions_cnt() + y_regions_cnt() - 1; ++diag){
        index_t x_begin, y_begin, nbr_elts;
        if(diag < x_regions_cnt()){ // Begin on the upper side
            x_begin = diag;
            y_begin = 0;
            nbr_elts = std::min(diag, y_regions_cnt());
        }
        else{ // Begin on the right side
            x_begin = x_regions_cnt() - 1;
            y_begin = diag - x_regions_cnt() + 1;
            nbr_elts = std::min(x_regions_cnt(), x_regions_cnt() + y_regions_cnt() - diag - 1);
        }
        index_t offs = 0;
        for(index_t nxt_offs = 1; nxt_offs < nbr_elts; ++nxt_offs){
            if(get_region(x_begin-nxt_offs, y_begin+nxt_offs).capacity() > 0){
                redo_bipartition(get_region(x_begin-offs, y_begin+offs), get_region(x_begin-nxt_offs, y_begin+nxt_offs));
                offs = nxt_offs;
            }
        }
    }

    
    // Perform optimization on the columns
    for(index_t x = 0; x < x_regions_cnt(); ++x){
        index_t y=0;
        for(index_t nxt_y = 1; nxt_y < y_regions_cnt(); ++nxt_y){
            if(get_region(x, nxt_y).capacity() > 0){
                redo_bipartition(get_region(x, y), get_region(x, nxt_y));
                y = nxt_y;
            }
        }
    }

    // Perform optimization on the rows
    for(index_t y = 0; y < y_regions_cnt(); ++y){
        index_t x=0;
        for(index_t nxt_x = 1; nxt_x < x_regions_cnt(); ++nxt_x){
            if(get_region(nxt_x, y).capacity() > 0){
                redo_bipartition(get_region(x, y), get_region(nxt_x, y));
                x = nxt_x;
            }
        }
    }
}

void region_distribution::region::redo_partition(std::vector<std::reference_wrapper<region_distribution::region> > regions){
    std::vector<cell_ref> all_cells;
    std::vector<capacity_t> caps;
    for(region_distribution::region & R : regions){
        all_cells.insert(all_cells.end(), R.cell_references_.begin(), R.cell_references_.end());
        caps.push_back(R.capacity_);
        R.cell_references_.clear();
    }
    std::sort(all_cells.begin(), all_cells.end(), [](cell_ref const a, cell_ref const b){ return a.allocated_capacity_ > b.allocated_capacity_; });

    current_allocation transporter(caps);
    for(auto const C : all_cells){
        std::vector<float_t> costs;
        for(region_distribution::region & R : regions){
            costs.push_back(R.distance(C) * static_cast<float_t>(C.allocated_capacity_));
        }
        transporter.add_source(C.allocated_capacity_, costs);
    }
    auto res = transporter.get_allocations();
    assert(res.size() == regions.size());
    for(index_t i=0; i<regions.size(); ++i){
        assert(res[i].size() == all_cells.size());
        for(index_t j=0; j<all_cells.size(); ++j){
            if(res[i][j] > 0){
                cell_ref C = all_cells[j];
                C.allocated_capacity_ = res[i][j];
                regions[i].get().cell_references_.push_back(C);
            }
        }
    }
}

void region_distribution::redo_quadpartitions(){
    for(index_t x = 1; x+1 < x_regions_cnt(); x+=2){
        for(index_t y = 1; y+1 < y_regions_cnt(); y+=2){
            region::redo_partition(std::vector<std::reference_wrapper<region> >({get_region(x,y), get_region(x,y+1), get_region(x+1,y), get_region(x+1,y+1)}));
        }
    }
    for(index_t x = 0; x+1 < x_regions_cnt(); x+=2){
        for(index_t y = 0; y+1 < y_regions_cnt(); y+=2){
            region::redo_partition(std::vector<std::reference_wrapper<region> >({get_region(x,y), get_region(x,y+1), get_region(x+1,y), get_region(x+1,y+1)}));
        }
    }
}

void region_distribution::fractions_minimization(){
    for(region & R : placement_regions_){
        R.uniquify_references();
    }

    // Find cycles of cut cells, then find a spanning tree to reallocate the cells
    // TODO
}

region_distribution::region_distribution(box<int_t> placement_area, std::vector<movable_cell> all_cells, std::vector<fixed_cell> all_obstacles) : x_cuts_cnt_(0), y_cuts_cnt_(0), placement_area_(placement_area), cell_list_(all_cells){

    std::vector<cell_ref> references;
    for(index_t i=0; i<all_cells.size(); ++i){
        movable_cell const & c = all_cells[i];
        if(c.demand_ == 0){
            throw std::runtime_error("A cell has been found with demand 0");
        }
        references.push_back( cell_ref(c.demand_, c.pos_, i) );
    }
    placement_regions_.push_back(
        region(placement_area_, all_obstacles, references)
    );
}

std::vector<region_distribution::movable_cell> region_distribution::export_positions() const{
    std::vector<point<float_t> > weighted_pos(cell_list_.size(), point<float_t>(0.0, 0.0));

    for(region const & R : placement_regions_){
        for(cell_ref C : R.cell_references_){
            weighted_pos[C.index_in_list_] = weighted_pos[C.index_in_list_] + static_cast<float_t>(C.allocated_capacity_) * R.pos_;
        }
    }

    std::vector<movable_cell> ret;
    for(index_t i=0; i<cell_list_.size(); ++i){
        movable_cell C = cell_list_[i];
        C.pos_ = ( static_cast<float_t>(1.0) / static_cast<float_t>(C.demand_) ) * weighted_pos[i];
        ret.push_back(C);
    }
    return ret;
}

struct OSRP_task{
    float_t size;
    float_t goal_pos;
    float_t weight;
    index_t orig;
    OSRP_task(float_t p, float_t s, float_t w, index_t i) : size(s), goal_pos(p), weight(w), orig(i) {}
    OSRP_task(){}
};

std::vector<float_t> get_optimal_quadratic_pos(std::vector<OSRP_task> cells, float_t pos_begin, float_t pos_end){

    if(cells.empty()){ return std::vector<float_t>(); }

    struct OSRP_cluster{
        index_t begin, end;
        float_t pos, size;
        float_t weight;
        OSRP_cluster(index_t b, index_t e, float_t p, float_t s, float_t w) : begin(b), end(e), pos(p), size(s), weight(w) {}
        void merge(OSRP_cluster const o){
            begin = o.begin;
            pos = (weight * pos + o.weight * o.pos) / (weight + o.weight);
            size += o.size;
            weight += o.weight;
        }
    };

    std::sort(cells.begin(), cells.end(), [](OSRP_task a, OSRP_task b){ return a.goal_pos < b.goal_pos; });

    // Modify the goal pos to get absolute goal positions between pos_begin and pos_end - sum_sizes
    float_t sum_sizes = 0.0;
    for(auto & c : cells){
        c.goal_pos -= sum_sizes;
        sum_sizes += c.size;
    }
    float_t abs_begin = pos_begin + 0.5 * cells[0].size; // First cell must be far enough from the beginning
    float_t abs_end = pos_end - sum_sizes + 0.5 * cells[0].size; // Last cell must be far enough from the end
    for(auto & c : cells){
        c.goal_pos = std::max(c.goal_pos, abs_begin);
        c.goal_pos = std::min(c.goal_pos, abs_end);
    }

    std::vector<OSRP_cluster> clusters;
    for(index_t i=0; i<cells.size(); ++i){
        OSRP_cluster to_add(
            i, i,
            cells[i].goal_pos,
            cells[i].size, cells[i].weight
        );
        while(not clusters.empty() and (clusters.back().pos >= to_add.pos)){
            to_add.merge(clusters.back());
            clusters.pop_back();
        }

        clusters.push_back(to_add);
    }

    std::vector<float_t> ret(cells.size(), 0.0);
    // For every cell, recover its true position from the absolute position of a cluster
    float_t sum_prev_sizes = 0.0;
    for(OSRP_cluster cur : clusters){
        for(index_t i=cur.begin; i <= cur.end; ++i){
            ret[cells[i].orig] = cur.pos + sum_prev_sizes;
            sum_prev_sizes += cells[i].size;
        }
    }
    return ret;
}




std::vector<region_distribution::movable_cell> region_distribution::export_spread_positions() const{
    std::vector<point<float_t> > weighted_pos(cell_list_.size(), point<float_t>(0.0, 0.0));

    for(region const & R : placement_regions_){
        index_t n = R.cell_references_.size();
        float_t total_capacity = static_cast<float_t>(R.capacity());
        box<float_t> surface = static_cast<box<float_t> >(R.surface_);

        std::vector<OSRP_task> x_cells(n), y_cells(n);
        for(index_t i=0; i<n; ++i){
            point<float_t> pt = R.cell_references_[i].pos_;
            float_t cap = static_cast<float_t>(R.cell_references_[i].allocated_capacity_);
            x_cells[i] = OSRP_task(pt.x_, cap/total_capacity * (surface.x_max_ - surface.x_min_), 1.0, i);
            y_cells[i] = OSRP_task(pt.y_, cap/total_capacity * (surface.y_max_ - surface.y_min_), 1.0, i);
        }
        std::vector<float_t> x_ret = get_optimal_quadratic_pos(x_cells, surface.x_min_, surface.x_max_);
        std::vector<float_t> y_ret = get_optimal_quadratic_pos(y_cells, surface.y_min_, surface.y_max_);

        for(index_t i=0; i<n; ++i){
            weighted_pos[R.cell_references_[i].index_in_list_] +=
                  static_cast<float_t>(R.cell_references_[i].allocated_capacity_)
                * point<float_t>(x_ret[i], y_ret[i]);
        }
    }

    std::vector<movable_cell> ret;
    for(index_t i=0; i<cell_list_.size(); ++i){
        movable_cell C = cell_list_[i];
        C.pos_ = ( static_cast<float_t>(1.0) / static_cast<float_t>(C.demand_) ) * weighted_pos[i];
        ret.push_back(C);
    }
    return ret;
}


} // Namespace gp
} // Namespace coloquinte

