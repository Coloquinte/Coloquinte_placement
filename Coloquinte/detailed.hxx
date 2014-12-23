
#ifndef COLOQUINTE_DETAILED
#define COLOQUINTE_DETAILED

#include "common.hxx"
#include "netlist.hxx"

#include <vector>
#include <limits>

namespace coloquinte{
namespace dp{

const index_t null_ind = std::numeric_limits<index_t>::max();

struct detailed_placement{
    struct internal_cell{
        point<int_t> position;
        int_t width;
        index_t row;

        // In order to get the neighbours in the detailed placement
        index_t neighbours_begin;
        index_t height;
        bool x_orientation, y_orientation;
    };

    // The placement region
    int_t min_x_, max_x_;
    int_t y_origin_;
    int_t row_height_;

    // Encode the topological state of the circuit: which cells are near each other
    // Makes extracting part of the circuit or optimizing positions at fixed topology easy
    std::vector<std::pair<index_t, index_t> > neighbours_; // The cells before and after on each row; cells spanning multiple columns use several positions

    std::vector<index_t> row_first_cells_, row_last_cells_; // For each row, which cells are the on the boundaries

    // All the cells and their properties
    std::vector<internal_cell> cells_;

    // Tests the coherency between positions, widths and topological representation
    void selfcheck() const;

    detailed_placement(
            std::vector<internal_cell> const cells,
            std::vector<std::vector<index_t> > const rows,
            int_t min_x, int_t max_x,
            int_t y_origin,
            index_t nbr_rows, int_t row_height
        );

    index_t cell_hght(index_t c) const{ return cells_[c].height; }
    index_t cell_cnt() const{ return cells_.size(); }
    index_t row_cnt()  const{ return row_first_cells_.size(); }
    index_t neighbour_index(index_t c, index_t r) const{
        assert( r >= cells_[c].row and r < cells_[c].row + cells_[c].height);
        return cells_[c].neighbours_begin + r - cells_[c].row;
    }
};

void optimize_positions(netlist const & circuit, detailed_placement & pl);
void optimize_swaps(netlist const & circuit, detailed_placement & pl, index_t row_extent, index_t cell_extent);

} // namespace dp
} // namespace coloquinte

#endif

