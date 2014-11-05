
#include "gp/solvers.hxx"

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>

namespace coloquinte{
namespace gp{

    std::vector<float_t> linear_system::solve_CG(std::vector<float_t> guess, float_t tol){
        Eigen::SparseMatrix<float_t> mat = Eigen::SparseMatrix<float_t>(size(), size());
        Eigen::ConjugateGradient<Eigen::SparseMatrix<float_t> > cg;

        std::vector<Eigen::Triplet<float_t> > triplets;
        for(auto t : matrix_){
            triplets.push_back(Eigen::Triplet<float_t>(t.r_, t.c_, t.val_));
        }
        mat.setFromTriplets(triplets.begin(), triplets.end());

        Eigen::Matrix<float_t, Eigen::Dynamic, 1> goal(size()), eguess(size());
        for(index_t i=0; i<size(); ++i){
            goal[i] = target_[i];
            eguess[i] = guess[i];
        }

        cg.setTolerance(tol / 1000000.0);
        cg.compute(mat);
        Eigen::Matrix<float_t, Eigen::Dynamic, 1> res(size());
        res = cg.solveWithGuess(goal, eguess);

        std::vector<float_t> ret;
        ret.resize(size());
        
        for(index_t i=0; i<size(); ++i){
            ret[i] = res[i];
        }
        
        return ret;
    }
/*
    std::vector<float_t> linear_system::solve_LLT();
    std::vector<float_t> linear_system::solve_LDLT();
*/

}
}



