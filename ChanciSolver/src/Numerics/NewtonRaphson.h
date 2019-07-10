#ifndef JACOBIAN_H
#define JACOBIAN_H

#include <Eigen/Sparse>
#include <vector>

template<typename FlowFunction, typename AccumulationFunction>
class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat;
    typedef Eigen::Triplet<double> Triplet;

    Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double, int> > solver;
    
    SparseMat jacobian;
    Eigen::VectorXd residual;
    Eigen::VectorXd solutiondelta;
    
 public:
    void iterate(FlowFunction calculateFlow, AccumulationFunction calculateAccumulation);
    void calculateResidual(FlowFunction calculateFlow, AccumulationFunction calculateAccumulation);
    void solve();
    void updateInteration(std::vector<double*> unknowns);
};

template<typename FlowFunction, typename AccumulationFunction>
    void NewtonRaphson<FlowFunction, AccumulationFunction>::solve(){
    solver.compute(jacobian);
    solutiondelta = solver.solve(residual);
};

template<typename FlowFunction, typename AccumulationFunction>
void NewtonRaphson<FlowFunction, AccumulationFunction>::updateInteration(std::vector<double*> unknowns){
    
};
template<typename FlowFunction, typename AccumulationFunction>
void NewtonRaphson<FlowFunction,AccumulationFunction>::iterate(FlowFunction calculateFlow, AccumulationFunction calculateAccumulation){

};

#endif /* JACOBIAN_H */
