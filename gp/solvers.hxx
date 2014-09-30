
#ifndef COLOQUINTE_GP_SOLVERS
#define COLOQUINTE_GP_SOLVERS

#include <vector>

namespace coloquinte{
namespace gp{

class linear_system{
    index_t mat_size; // Results up to mat_size are used by the placement
    index_t private_mat_size; // The vector and matrix have up to private_mat_size rows, but the rows after mat_size are only intermediate variables

    struct mat_elt{
        index_t x_ind, y_ind;
        float_t weight;

        mat_elt(index_t x_index, index_t y_index, float_t delta_weight) : x_ind(x_index), y_ind(y_index), weight(delta_weight){}
    };
    struct vec_elt{
        index_t ind;
        float_t weight;

        vec_elt(index_t y_index, float_t delta_weight) : ind(y_index), weight(delta_weight){}
    };

    std::vector<mat_elt> Mat;
    std::vector<vec_elt> Vec;
    
    public:
    void add_to_mat(index_t x_index, index_t y_index, float_t delta_weight){ Mat.push_back(mat_elt(x_index, y_index, delta_weight)); }
    void add_to_vec(index_t y_index, float_t delta_weight){ Vec.push_back(vec_elt(y_index, delta_weight)); }
};

class csr_sys{
    struct csr_elt{
        index_t x_ind;
        float_t value;
    };

    std::vector<index_t> line_begins;
    std::vector<csr_elt> elements;

    std::vector<float_t> Vec;

    std::vector<float_t> cg(FLOAR err);
};

} // namespace gp
} // namespace coloquinte

#endif


