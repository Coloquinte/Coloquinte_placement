
#include "gp/rough_legalizers.hxx"
#include <algorithm>
#include <cmath>
#include <cassert>

namespace coloquinte{
namespace gp{

void region_distribution::region::selfcheck() const{
    capacity_t total_allocated = 0;
    for(cell_ref const c : cell_references_){
        total_allocated += c.allocated_capacity_;
        assert(c.allocated_capacity_ > 0);
    }
    assert(total_allocated + unused_capacity_ == capacity_);
}

void region_distribution::selfcheck() const{
    for(region const & R : placement_regions_){
        R.selfcheck();
    }
    std::vector<capacity_t> capacities(cell_list_.size(), 0);
    for(region const & R : placement_regions_){
        for(cell_ref const C : R.cell_references_){
            capacities[C.index_in_list_] += C.allocated_capacity_;
        }
    }
    for(index_t i=0; i < cell_list_.size(); ++i){
        assert(capacities[i] == cell_list_[i].demand_);
    }
}

region_distribution::region::region(box<int_t> bx, std::vector<fixed_cell> obstacles, std::vector<cell_ref> cells) : surface_(bx), cell_references_(cells), obstacles_(obstacles){
    x_pos_ = static_cast<float_t>(surface_.x_max_ + surface_.x_min_) * 0.5;
    y_pos_ = static_cast<float_t>(surface_.y_max_ + surface_.y_min_) * 0.5;

    capacity_ = static_cast<capacity_t>(surface_.x_max_ - surface_.x_min_) * static_cast<capacity_t>(surface_.y_max_ - surface_.y_min_);
    for(auto const & O : obstacles_){
        // TODO: do it properly
    }

    unused_capacity_ = capacity_;
    for(auto const & C : cell_references_){
        assert(unused_capacity_ >= C.allocated_capacity_);
        unused_capacity_ -= C.allocated_capacity_;
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
    assert(lft.unused_capacity() + rgt.unused_capacity() == unused_capacity());
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
    assert(up.unused_capacity() + dwn.unused_capacity() == unused_capacity());
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

    std::vector<cell_ref> cells_a_side, cells_b_side;
    if(preference_limit >= b_capacity_limit and preference_limit <= a_capacity_limit){
        cells_a_side.insert(cells_a_side.end(), cells.begin(), cells.begin() + preference_limit);
        cells_b_side.insert(cells_b_side.end(), cells.begin() + preference_limit, cells.end());
    }
    else{
        index_t cut_position;
        capacity_t allocated_to_a_part;
        if(preference_limit < b_capacity_limit){ // Pack on b
            cut_position = b_capacity_limit-1; // Exists since preference_limit >= 0
            allocated_to_a_part = cells[cut_position].allocated_capacity_ - remaining_capacity_b;
        }
        else{ // Pack on a
        // if(preference_limit > a_capacity_limit)
            cut_position = a_capacity_limit; // Necessarily a correct position since preference_limits <= cells.size()
            allocated_to_a_part = remaining_capacity_a;
        }

        cells_a_side.insert(cells_a_side.end(), cells.begin(), cells.begin() + cut_position);

        cell_ref cell_cut_a = cells[cut_position], cell_cut_b = cells[cut_position];
        cell_cut_a.allocated_capacity_ = allocated_to_a_part;
        cell_cut_b.allocated_capacity_ -= allocated_to_a_part;
        if(cell_cut_a.allocated_capacity_ > 0){ cells_a_side.push_back(cell_cut_a); }
        if(cell_cut_b.allocated_capacity_ > 0){ cells_b_side.push_back(cell_cut_b); }

        cells_b_side.insert(cells_b_side.end(), cells.begin() + cut_position+1, cells.end());
    }

    capacity_t unused_capacity_a = region_a.capacity_;
    capacity_t unused_capacity_b = region_b.capacity_;

    for(cell_ref const & c : cells_a_side){
        unused_capacity_a -= c.allocated_capacity_;
    }
    for(cell_ref const & c : cells_b_side){
        unused_capacity_b -= c.allocated_capacity_;
    }

    
    // Verify that the cells have been correctly handled
    std::vector<float_t> costs_a_side, costs_b_side;
    for(auto const C : cells_a_side){ costs_a_side.push_back(C.marginal_cost_); }
    for(auto const C : cells_b_side){ costs_b_side.push_back(C.marginal_cost_); }
    if(not costs_a_side.empty() and not costs_b_side.empty()){
        float_t max_left = *std::max_element(costs_a_side.begin(), costs_a_side.end()), min_right = *std::min_element(costs_b_side.begin(), costs_b_side.end());
        assert(max_left <= min_right);
    }
    

    region_a.cell_references_ = cells_a_side;
    region_b.cell_references_ = cells_b_side;

    region_a.unused_capacity_ = unused_capacity_a;
    region_b.unused_capacity_ = unused_capacity_b;

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
    selfcheck();
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
    selfcheck();
}


void region_distribution::quadpartition(){
    // Naive implementation; maybe advanced ones could use an exact solution
    x_bipartition();
    y_bipartition();
}

float_t region_distribution::cost() const{
    double res = 0.0;
    for(region const & R : placement_regions_){
        double this_region = 0.0;
        for(cell_ref const & C : R.cell_references_){
            this_region += static_cast<double>(C.allocated_capacity_) * static_cast<double>(R.distance(C));
        }
        res += this_region;
    }
    return res;
}

void region_distribution::redo_bipartition(region_distribution::region & Ra, region_distribution::region & Rb){
    std::vector<cell_ref> cells;
    for(cell_ref C : Ra.cell_references_){ cells.push_back(C); }
    for(cell_ref C : Rb.cell_references_){ cells.push_back(C); }

    region::distribute_new_cells(Ra, Rb, cells);
}

void region_distribution::redo_bipartitions(){
    for(index_t x = 0; x+1 < x_regions_cnt(); ++x){
        for(index_t y = 0; y+1 < y_regions_cnt(); ++y){
            redo_bipartition(get_region(x,y), get_region(x+1, y));
            redo_bipartition(get_region(x,y), get_region(x, y+1));
            redo_bipartition(get_region(x,y), get_region(x+1, y+1));
        }
    }
    selfcheck();
}

/*
region_distribution::region_distribution(global_placement const & orig) : x_cuts_cnt_(0), y_cuts_cnt_(0){
    placement_regions_.resize(1);
    
    std::vector<movable_cell> cells;
    std::vector<cell_ref> cell_refs;
    std::vector<fixed_cell> obstacles;
    for(index_t i=0; i<orig.placement_.size(); ++i){

        if(orig.reference_circuit_->get_cell(i).movability_.is_fully_movable()){
            movable_cell c;
            c.demand_ = orig.reference_circuit_->get_cell(i).area_;
            c.x_pos_ = orig.placement_[i].x_pos_;
            c.y_pos_ = orig.placement_[i].y_pos_;
            c.index_in_placement_ = i;

            cell_ref cr;
            cr.allocated_capacity_ = orig.reference_circuit_->get_cell(i).area_;
            cr.x_pos_ = orig.placement_[i].x_pos_;
            cr.y_pos_ = orig.placement_[i].y_pos_;
            cr.index_in_list_ = i;

            cells.push_back(c);
            cell_refs.push_back(cr);
        }
        else{
            float_t x_min = orig.placement_[i].x_pos_ - orig.reference_circuit_->get_cell(i).x_size_ * 0.5,
                    x_max = orig.placement_[i].x_pos_ + orig.reference_circuit_->get_cell(i).x_size_ * 0.5,
                    y_min = orig.placement_[i].y_pos_ - orig.reference_circuit_->get_cell(i).y_size_ * 0.5,
                    y_max = orig.placement_[i].y_pos_ + orig.reference_circuit_->get_cell(i).y_size_ * 0.5;
            obstacles.push_back(fixed_cell(box<int_t>(x_min, x_max, y_min, y_max)));
        }
    }

    placement_area_ = orig.reference_circuit_->placement_area();
    placement_regions_[0] = region(placement_area_, obstacles, cell_refs);
}
*/

region_distribution::region_distribution(box<int_t> placement_area, std::vector<movable_cell> all_cells) : placement_area_(placement_area), x_cuts_cnt_(0), y_cuts_cnt_(0), cell_list_(all_cells){

    std::vector<cell_ref> references;
    for(index_t i=0; i<all_cells.size(); ++i){
        movable_cell const & c = all_cells[i];
        references.push_back( cell_ref(c.demand_, c.x_pos_, c.y_pos_, i) );
    }

    placement_regions_.push_back(
        region(placement_area_, std::vector<fixed_cell>(), references)
    );

    selfcheck();
}



} // Namespace gp
} // Namespace coloquinte

