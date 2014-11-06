
#ifndef COLOQUINE_GP_SOLVERS
#define COLOQUINE_GP_SOLVERS

#include "common.hxx"

#include <vector>

namespace coloquinte{
namespace gp{

struct matrix_triplet{
    index_t r_, c_;
    float_t val_;
    matrix_triplet(index_t ri, index_t ci, float_t v) : r_(ri), c_(ci), val_(v){}
    bool operator<(matrix_triplet const o){ return r_ < o.r_ || (r_ == o.r_ && c_ < o.c_); }
};

class linear_system{
    std::vector<matrix_triplet> matrix_;
    std::vector<float_t> target_;

    public:
    void add_triplet(index_t row, index_t col, float_t val){ matrix_.push_back(matrix_triplet(row, col, val)); }

    linear_system operator+(linear_system const & o){
        if(o.target_.size() != target_.size()){ throw std::runtime_error("Mismatched system size"); }
        linear_system ret(size());

        ret.matrix_ = matrix_;
        ret.matrix_.insert(ret.matrix_.end(), o.matrix_.begin(), o.matrix_.end());

        ret.target_.resize(target_.size());
        for(index_t i=0; i<ret.target_.size(); ++i){
            ret.target_[i] = target_[i] + o.target_[i];
        }

        return ret;
    }

    void add_doublet(index_t row, float_t val){
        target_[row] += val;
    }

    void add_force(
        float_t tol,
        float_t scale,
        index_t c1,    index_t c2,
        float_t pos1,  float_t pos2,
        float_t offs1, float_t offs2
    ){
        float_t force = scale / std::max(std::abs(pos1 - pos2), tol);
        add_triplet(c1, c1, force);
        add_triplet(c2, c2, force);
        add_triplet(c1, c2, -force);
        add_triplet(c2, c1, -force);
        add_doublet(c1, force * (offs2-offs1));
        add_doublet(c2, force * (offs1-offs2));
    }

    void add_fixed_force(
        float_t tol,
        float_t scale,
        index_t c,
        float_t cell_pos, float_t fixed_pos,
        float_t offs
    ){
        float_t force = scale / std::max(std::abs(cell_pos - fixed_pos), tol);
        add_triplet(c, c, force);
        add_doublet(c, force * (fixed_pos-offs));
    }

    linear_system(index_t s) : target_(s, 0.0){}

    index_t size() const{ return target_.size(); }

    std::vector<float_t> solve_CG(std::vector<float_t> guess, float_t tol);
    std::vector<float_t> solve_cholesky();
    

};

} // namespace gp
} // namespace coloquinte

#endif


