
#include "Coloquinte/detailed.hxx"

#include <cassert>

namespace coloquinte{
namespace dp{

detailed_placement::detailed_placement(
        std::vector<int_t>    positions,
        std::vector<index_t>  row_positions,
        std::vector<int_t>    widths,
        std::vector<index_t>  heights,
        std::vector<bool> x_orientations,
        std::vector<bool> y_orientations,
        std::vector<std::vector<index_t> > const rows,
        int_t min_x, int_t max_x,
        int_t y_origin,
        index_t nbr_rows, int_t row_height
    )
    :
        min_x_(min_x), max_x_(max_x),
        y_origin_(y_origin),
        row_height_(row_height),
        positions_(positions),
        cell_rows_(row_positions),
        widths_(widths),
        x_orientations_(x_orientations),
        y_orientations_(y_orientations)
    {

    index_t sz = positions.size();

    assert(row_height > 0);
    assert(min_x < max_x);
    assert(row_positions.size() == sz
        and widths.size()  == positions.size()
        and heights.size() == positions.size()
    );
    assert(rows.size() == nbr_rows);

    cell_lims_.resize(1, 0);
    index_t nbr_lims = 0; 
    for(index_t h : heights){
        nbr_lims += h;
        cell_lims_.push_back(nbr_lims);
    }

    cells_before_ .resize(nbr_lims, null_ind);
    cells_after_  .resize(nbr_lims, null_ind);

    row_first_cells_ .resize(nbr_rows, null_ind);
    row_last_cells_  .resize(nbr_rows, null_ind);

    std::vector<bool> explored(nbr_lims, false);
    // Now we extract the dependencies
    for(index_t r=0; r<rows.size(); ++r){

        if(not rows[r].empty()){
            row_first_cells_[r] = rows[r].front();
            row_last_cells_[r]  = rows[r].back();
        }

        for(index_t c : rows[r]){
            // The row_positions are correct
            assert( r >= cell_rows_[c] and r < cell_rows_[c] + cell_hght(c));

            // Has this row of the cell already been visited?
            assert(not explored[cell_lims_[c] + r - cell_rows_[c]]);
            explored[cell_lims_[c] + r - cell_rows_[c]] = true;
        }

        for(index_t i=0; i+1<rows[r].size(); ++i){
            index_t c1 = rows[r][i], c2 = rows[r][i+1];
            assert(c1 < 1000000 && c2 < 1000000);

            // Save in the internal format
            cells_after_ [cell_lims_[c1] + r - cell_rows_[c1]] = c2;
            cells_before_[cell_lims_[c2] + r - cell_rows_[c2]] = c1;

            // The positions are correct
            assert( positions[c1] + widths[c1] <= positions[c2] );
        }
    }

    // Every level of every cell must have been visited
    for(bool o : explored)
        assert(o);

    // Verify that we haven't made any obvious mistake
    selfcheck();
}

void detailed_placement::selfcheck() const{
    assert(cells_before_.size()   == cells_after_.size());
    assert( cell_rows_.size()     == cell_cnt()
        and widths_.size()        == cell_cnt()
        and positions_.size()     == cell_cnt()
        and x_orientations_.size() == cell_cnt()
        and y_orientations_.size() == cell_cnt()
    );
    assert(row_first_cells_.size() == row_last_cells_.size());

    for(index_t i=0; i<cell_cnt(); ++i){
        for(index_t l=cell_lims_[i]; l<cell_lims_[i+1]; ++l){
            // Correct x position
            assert(positions_[i] >= min_x_ and positions_[i] + widths_[i] <= max_x_);
            assert(cell_rows_[i] + cell_hght(i) <= row_cnt());

            if(cells_before_[l] != null_ind){
                // Correct neighbour position
                assert(positions_[i] >= positions_[cells_before_[l]] + widths_[cells_before_[l]]);
            }
            else{
                // Beginning of a row
                assert(row_first_cells_[cell_rows_[i] + l - cell_lims_[i]] == i);
            }
            if(cells_after_[l] != null_ind){
                // Correct neighbour position
                assert(positions_[i] + widths_[i] <= positions_[cells_after_[l]]);
            }
            else{
                // End of a row
                assert(row_last_cells_[cell_rows_[i] + l - cell_lims_[i]] == i);
            }
        }
    }
    
    for(index_t i=0; i<row_cnt(); ++i){
        if(row_first_cells_[i] != null_ind){
            assert(cells_before_[cell_lims_[row_first_cells_[i]] + i - cell_rows_[row_first_cells_[i]]] == null_ind);
        }
        if(row_last_cells_[i] != null_ind){
            assert(cells_after_[cell_lims_[row_last_cells_[i]] + i - cell_rows_[row_last_cells_[i]]] == null_ind);
        }
    }
}

} // namespace dp
} // namespace coloquinte

