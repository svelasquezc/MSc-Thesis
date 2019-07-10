#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <Eigen/Sparse>
#include <vector>

#include "Mesh.h"
#include "Rock.h"
#include "Equilibrium_Relation.h"

template<typename FlowFunction, typename AccumulationFunction>
class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat;
    typedef Eigen::Triplet<double> Triplet;
    
    Eigen::BiCGSTAB<SparseMat, Eigen::IncompleteLUT<double, int> > solver;
    
    SparseMat jacobian;
    Eigen::VectorXd residual;
    Eigen::VectorXd initial_residual;
    Eigen::VectorXd solutiondelta;

    FlowFunction calculateFlow;
    AccumulationFunction calculateAccumulation;

    int iteration=0;
    
 public:
    void NewtonRaphson(FlowFunction _calculateFlow, AccumulationFunction _calculateAccumulation) :
    calculateFlow(_calculateFlow), calculateAccumulation(_calculateAccumulation){};
    
    void iterate(const double relative_change_in_residual, const int& _term, const Mesh& _mesh, std::vector<std::shared_ptr<Fluid>>& _characterized_fluids, Rock& _rock);
    double calculateResidual(const int& _term, std::vector<std::shared_ptr<Fluid>>& _characterized_fluids, Mesh& _mesh, Rock& _rock);
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
void NewtonRaphson<FlowFunction,AccumulationFunction>::iterate(const double relative_change_in_residual, const int& _term, const Mesh& _mesh, std::vector<std::shared_ptr<Fluid>>& _characterized_fluids, Rock& _rock){
    
    do{

        for(auto fluid : _characterized_fluids){
            
        }
        
    }while(residual.squaredNorm()/initial_residual.squaredNorm() > relative_change_in_residual);
    
};

#endif /* NEWTONRAPHSON_H */
