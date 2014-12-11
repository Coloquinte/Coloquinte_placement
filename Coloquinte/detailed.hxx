
#ifndef COLOQUINTE_DETAILED
#define COLOQUINTE_DETAILED

#include "common.hxx"
#include "netlist.hxx"

#include <vector>

namespace coloquinte{
namespace dp{

const index_t null_ind = 0;

struct detailed_placement{
    // The placement region
    int_t min_x_, max_x_;
    int_t y_origin_;
    int_t row_height_;

    // Encode the topological state of the circuit: which cells are near each other
    // Makes extracting part of the circuit or optimizing positions at fixed topology easy
    std::vector<index_t> cell_lims_; // The limit difference gives the number of rows occupied by a cell
    std::vector<index_t> cells_before_, cells_after_; // The cells before and after on each row; cells spanning multiple columns use several positions

    std::vector<index_t> row_first_cells_, row_last_cells_; // For each row, which cells are the on the boundaries

    // The current position of a cell
    std::vector<int_t> positions_;
    // The current row of a cell (smallest row number occupied)
    std::vector<index_t> cell_rows_;
    // The size of the cell
    std::vector<int_t> widths_;

    std::vector<bool> x_orientations_, y_orientations_;
    // Tests the coherency between positions, widths and topological representation
    void selfcheck() const;

    detailed_placement(
            std::vector<int_t> positions,
            std::vector<index_t> row_positions,
            std::vector<int_t> widths,
            std::vector<index_t> heights,
            std::vector<bool> x_orientations,
            std::vector<bool> y_orientations,
            std::vector<std::vector<index_t> > rows,
            int_t min_x, int_t max_x,
            int_t y_origin,
            index_t nbr_rows, int_t row_height
        );

    index_t cell_hght(index_t c) const{ return cell_lims_[c+1] - cell_lims_[c]; }
    index_t cell_cnt() const{ return cell_lims_.size() - 1; }
    index_t row_cnt()  const{ return row_first_cells_.size(); }
};

float_t get_HPWL_wirelength(netlist const & circuit, detailed_placement const & pl);

} // namespace dp
} // namespace coloquinte

#endif

