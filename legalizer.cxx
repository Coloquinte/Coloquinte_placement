
#include "Coloquinte/legalizer.hxx"

#include <algorithm>
#include <cmath>

namespace coloquinte{
namespace dp{

void get_result(netlist const & circuit, detailed_placement const & dpl, gp::placement_t & gpl){
    for(index_t c=0; c<circuit.cell_cnt(); ++c){
        if( (circuit.get_cell(c).attributes & XMovable) != 0)
            gpl.positions_[c].x_ = static_cast<float_t>(dpl.positions_[c]) + 0.5 * circuit.get_cell(c).size.x_;
        if( (circuit.get_cell(c).attributes & YMovable) != 0)
            gpl.positions_[c].y_ = static_cast<float_t>(dpl.y_origin_) + dpl.cell_rows_[c] * dpl.row_height_ + 0.5 * circuit.get_cell(c).size.y_;

        if( (circuit.get_cell(c).attributes & XFlippable) != 0)
            gpl.orientations_[c].x_ = dpl.x_orientations_[c] ? 1.0 : -1.0;
        if( (circuit.get_cell(c).attributes & YFlippable) != 0)
            gpl.orientations_[c].y_ = dpl.y_orientations_[c] ? 1.0 : -1.0;
    }
}

struct cell_to_leg{
    int_t x_pos, y_pos;
    index_t original_cell;
    int_t width;
    index_t nbr_rows;

    bool operator<(cell_to_leg const o) const{ return x_pos < o.x_pos; }

    cell_to_leg(int_t x, int_t y, index_t ind, int_t w, index_t rows)
    : x_pos(x), y_pos(y),
    original_cell(ind),
    width(w),
    nbr_rows(rows)
    {}
};

struct fixed_cell_interval{
    int_t min_x, max_x;
    index_t cell_ind;

    bool operator<(fixed_cell_interval const o) const{ return min_x > o.min_x; }
    fixed_cell_interval(int_t mn, int_t mx, index_t ind) : min_x(mn), max_x(mx), cell_ind(ind){}
};

detailed_placement legalize(netlist const & circuit, gp::placement_t const & pl, box<int_t> surface, int_t row_height){
    if(row_height <= 0) throw std::runtime_error("The rows' height should be positive\n");

    index_t nbr_rows = (surface.y_max_ - surface.y_min_) / row_height;
    // The position of the ith row is surface.y_min_ + i * row_height

    std::vector<std::vector<fixed_cell_interval> > row_occupation(nbr_rows);
    std::vector<cell_to_leg> cells;

    std::vector<index_t> row_index(circuit.cell_cnt());
    std::vector<int_t>   cell_pos(circuit.cell_cnt());
    std::vector<int_t>   cell_width(circuit.cell_cnt());
    std::vector<index_t> cell_height(circuit.cell_cnt());

    for(index_t i=0; i<circuit.cell_cnt(); ++i){
        auto cur = circuit.get_cell(i);
        // Assumes fixed if not both XMovable and YMovable
        if( (cur.attributes & XMovable) != 0 && (cur.attributes & YMovable) != 0){
            // Just truncate the position we target
            point<int_t> target_pos = pl.positions_[i] - static_cast<float_t>(0.5) * static_cast<point<float_t> >(cur.size);
            index_t cur_cell_rows = (cur.size.y_ + row_height -1) / row_height;
            cells.push_back(cell_to_leg(target_pos.x_, target_pos.y_, i, cur.size.x_, cur_cell_rows));
            cell_height[i] = cur_cell_rows;
            cell_width[i] = cur.size.x_;
        }
        else{
            // In each row, we put the index of the fixed cell and the range that is already occupied
            int_t low_x_pos  = std::round(pl.positions_[i].x_ - 0.5 * static_cast<float_t>(cur.size.x_)),
                  hgh_x_pos  = std::round(pl.positions_[i].x_ + 0.5 * static_cast<float_t>(cur.size.x_)),
                  low_y_pos  = std::round(pl.positions_[i].y_ - 0.5 * static_cast<float_t>(cur.size.y_)),
                  hgh_y_pos  = std::round(pl.positions_[i].y_ + 0.5 * static_cast<float_t>(cur.size.y_));

            assert(low_x_pos < hgh_x_pos and low_y_pos < hgh_y_pos);
            if(hgh_y_pos <= surface.y_min_ or low_y_pos >= surface.y_max_ or hgh_x_pos <= surface.x_min_ or low_x_pos >= surface.x_max_) continue; // No intersection

            hgh_x_pos = std::min(surface.x_max_, hgh_x_pos);
            hgh_y_pos = std::min(surface.y_max_, hgh_y_pos);
            low_x_pos = std::max(surface.x_min_, low_x_pos);
            low_y_pos = std::max(surface.y_min_, low_y_pos);

            index_t first_row = (low_y_pos - surface.y_min_) / row_height;
            index_t last_row = (index_t) (hgh_y_pos - surface.y_min_ + row_height - 1) / row_height; // Exclusive: if the cell spans the next row, i.e. pos % row_height >= 0, include it too
            assert(last_row <= nbr_rows);

            row_index[i] = first_row;
            cell_pos[i] = low_x_pos;
            cell_height[i] = last_row - first_row;
            cell_width[i] = hgh_x_pos - low_x_pos; // Only the part within the placement region
            for(index_t r=first_row; r<last_row; ++r){
                row_occupation[r].push_back(fixed_cell_interval(low_x_pos, hgh_x_pos, i));
            }
        }
    }

    for(std::vector<fixed_cell_interval> & L : row_occupation){
        std::sort(L.begin(), L.end()); // Sorts from last to first, so that we may use pop_back()
        // Doesn't collapse them yet, which may make for bigger complexities
        for(index_t i=0; i+1<L.size(); ++i){
            if(L[i].min_x < L[i+1].max_x)
                throw std::runtime_error("Sorry, I don't handle overlapping fixed cells yet\n");
        }
    }

    std::vector<std::vector<index_t> > cells_by_rows(nbr_rows);
    std::vector<int_t> first_available_position(nbr_rows, surface.x_min_);

    // Sort the cells by x position
    std::sort(cells.begin(), cells.end());

    for(cell_to_leg C : cells){
        // Dumb, quick and dirty best-fit legalization
        bool found_location = false;

        // Properties of the current best solution
        int_t best_x=0;
        int_t best_cost=0;
        index_t best_row=0;

        // Helper function
        auto check_row_cost = [&](index_t r, cell_to_leg const cell, int_t additional_cost){
            // Find where to put the cell in these rows
            // Simple method: get a range where we can put the cell

            assert(r + cell.nbr_rows <= nbr_rows);
            assert(additional_cost >= 0);

            // First position where we can put it
            int_t cur_pos = *std::max_element(first_available_position.begin() + r, first_available_position.begin() + r + cell.nbr_rows);
            int_t max_lim = surface.x_max_ - cell.width;
            int_t interval_lim;
            do{
                interval_lim = max_lim;
                // For each row, test if obstacles prevent us from putting a cell here
                // Until we find a correct position or are beyond the maximum position
                for(index_t i = 0; i<cell.nbr_rows; ++i){
                    // Find the first obstacle which is after this position
                    // TODO: use lower/upper bound
                    auto it=row_occupation[r+i].rbegin();
                    for(; it != row_occupation[r+i].rend() && it->max_x <= cur_pos; ++it){
                    }
                    if(it != row_occupation[r+i].rend()){ // There is an obstacle on the right
                        assert(it->min_x < it->max_x);
                        int_t cur_lim = it->min_x - cell.width; // Where the obstacles contrains us
                        interval_lim = std::min(cur_lim, interval_lim); // Constraint
                        if(cur_lim < cur_pos){ // If this particular obstacle constrained us so that it is not possible to make it here, we increment the position
                            cur_pos = std::max(it->max_x, cur_pos);
                        }
                    }
                }
                // Do it again until we are past the optimal solution
            }while(interval_lim < cur_pos and interval_lim < max_lim and cur_pos < max_lim); // Not admissible and we encountered an obstacle and there is still hope

            if(interval_lim >= cur_pos){ // An admissible solution is found
                // TODO: if the solution may not be the optimal one, try it again
                int_t row_best_x = std::min(interval_lim, std::max(cur_pos, cell.x_pos));
                int_t row_cost_x = std::abs(row_best_x - cell.x_pos);
                if(not found_location or row_cost_x + additional_cost < best_cost){
                    found_location = true;
                    best_cost = row_cost_x + additional_cost;
                    best_x = row_best_x;
                    best_row = r;
                }
            }
        };

        // The row where we would prefer the cell to go
        if(C.nbr_rows > nbr_rows) throw std::runtime_error("Impossible to legalize a cell spanning more rows than are available\n");
        index_t central_row = std::min( (index_t) std::max( (C.y_pos - surface.y_min_) / row_height, 0), nbr_rows-C.nbr_rows);

        // Try every possible row from the center, until we can't improve the cost
        for(index_t row_dist = 0;
            (central_row + row_dist < nbr_rows or central_row >= row_dist)
            and (not found_location or static_cast<int_t>(row_dist) * row_height < best_cost);
            ++row_dist
        ){
            if(central_row + row_dist < nbr_rows - C.nbr_rows){
                check_row_cost(central_row + row_dist, C, row_dist * row_height);
            }
            if(central_row >= row_dist){
                check_row_cost(central_row - row_dist, C, row_dist * row_height);
            }
        }

        if(not found_location){
            throw std::runtime_error("Didn't manage to pack a cell due to dumb algorithm\n");
        }
        else{
            assert(best_x + C.width <= surface.x_max_ and best_x >= surface.x_min_);
            // Update the occupied rows
            for(index_t r = best_row; r < best_row + C.nbr_rows; ++r){
                // Include the obstacles
                while(not row_occupation[r].empty()
                    and row_occupation[r].back().max_x <= best_x){
                    cells_by_rows[r].push_back(row_occupation[r].back().cell_ind);
                    row_occupation[r].pop_back();
                }
                assert(row_occupation[r].empty() or row_occupation[r].back().min_x >= best_x + C.width);

                cells_by_rows[r].push_back(C.original_cell);
                first_available_position[r] = best_x + C.width;
            }
            cell_pos[C.original_cell] = best_x;
            row_index[C.original_cell] = best_row;
        }
    }

    // Finally, push the remaining fixed cells
    for(index_t r=0; r<nbr_rows; ++r){
        while(not row_occupation[r].empty()){
            cells_by_rows[r].push_back(row_occupation[r].back().cell_ind);
            row_occupation[r].pop_back();
        }
    }

    std::vector<bool> x_orientation(circuit.cell_cnt()), y_orientation(circuit.cell_cnt());
    for(index_t c=0; c<circuit.cell_cnt(); ++c){
        x_orientation[c] = pl.orientations_[c].x_ >= 0.0;
        y_orientation[c] = pl.orientations_[c].y_ >= 0.0;
    }

    return detailed_placement(
        cell_pos,
        row_index,
        cell_width,
        cell_height,
        x_orientation,
        y_orientation,
        cells_by_rows,
        surface.x_min_, surface.x_max_,
        surface.y_min_,
        nbr_rows, row_height
    );
}

} // namespace dp
} // namespace coloquinte

