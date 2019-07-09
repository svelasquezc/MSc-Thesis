#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <Eigen/Sparse>
#include <vector>

class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat;
    typedef Eigen::Triplet<double> Triplet;

    Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double, int> > solver;
    
    SparseMat jacobian;
    Eigen::VectorXd residual;
    Eigen::VectorXd solutiondelta;
    
 public:
    template<typename FlowFunction, typename AccumulationFunction>
    void iterate(FlowFunction calculateFlow, AccumulationFunction calculateAccumulation);
    void solve();
    void updateInteration(std::vector<double*> unknowns);
};

void NewtonRaphson::solve(){
    solver.compute(jacobian);
    solutiondelta = solver.solve(residual);
};

void NewtonRaphson::updateInteration(std::vector<double*> unknowns){
    
};

template<typename FlowFunction, typename AccumulationFunction>
    void NewtonRaphson::iterate(FlowFunction calculateFlow, AccumulationFunction calculateAccumulation){

};

#endif /* JACOBIAN_H */
