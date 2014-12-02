
#include "gp/rough_legalizers.hxx"
#include "gp/transportation.hxx"

#include <algorithm>
#include <cmath>
#include <cassert>

namespace coloquinte{
namespace gp{

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
    x_regions_cnt_ *= 2;

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
    y_regions_cnt_ *= 2;

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
    redo_multipartitions(2);
}

void region_distribution::redo_multipartitions(index_t width){
    redo_multipartitions(width, width);
}

void region_distribution::redo_multipartitions(index_t x_width, index_t y_width){
    if(x_width < 2 or y_width < 2) throw std::runtime_error("Multipartitioning requires an optimization window of 2 or more\n");
    for(index_t x = x_width/2; x+x_width <= x_regions_cnt(); x+=x_width){
        for(index_t y = y_width/2; y+y_width <= y_regions_cnt(); y+=y_width){
            std::vector<std::reference_wrapper<region> > to_opt;
            for(index_t l_x=x; l_x < x + x_width; ++l_x){
                for(index_t l_y=y; l_y < y + y_width; ++l_y){
                    to_opt.push_back(std::reference_wrapper<region>(get_region(l_x, l_y)));
                }
            }
            region::redo_partition(to_opt);
        }
    }
    for(index_t x = 0; x+x_width <= x_regions_cnt(); x+=x_width){
        for(index_t y = 0; y+y_width <= y_regions_cnt(); y+=y_width){
            std::vector<std::reference_wrapper<region> > to_opt;
            for(index_t l_x=x; l_x < x + x_width; ++l_x){
                for(index_t l_y=y; l_y < y + y_width; ++l_y){
                    to_opt.push_back(std::reference_wrapper<region>(get_region(l_x, l_y)));
                }
            }
            region::redo_partition(to_opt);
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

region_distribution::region_distribution(box<int_t> placement_area, std::vector<movable_cell> all_cells, std::vector<fixed_cell> all_obstacles) : x_regions_cnt_(1), y_regions_cnt_(1), placement_area_(placement_area), cell_list_(all_cells){

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

