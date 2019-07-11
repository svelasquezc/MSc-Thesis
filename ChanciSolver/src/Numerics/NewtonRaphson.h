#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <Eigen/Sparse>
#include <vector>

#include "Mesh.h"
#include "Rock.h"
#include "Equilibrium_Relation.h"

template<typename FlowFunction_t, typename AccumulationFunction_t>
    class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat_t;
    typedef Eigen::Triplet<double> Triplet_t;
    
    Eigen::BiCGSTAB<SparseMat_t, Eigen::IncompleteLUT<double, int> > _solver;
    
    SparseMat_t _jacobian;
    Eigen::VectorXd _residual;
    Eigen::VectorXd _initial_residual;
    Eigen::VectorXd _solution_delta;

    FlowFunction_t _calculateFlow;
    AccumulationFunction_t _calculateAccumulation;

    int iteration=0;
    
 public:
    NewtonRaphson(FlowFunction_t calculateFlow, AccumulationFunction_t calculateAccumulation) :
    _calculateFlow(calculateFlow), _calculateAccumulation(calculateAccumulation){};
    
    void iterate(const double relative_change_in_residual, const int term, const Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, Rock& rock);
    double calculateResidual(const int& _term, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, Mesh& mesh, Rock& rock);
    void solve();
    void updateInteration(std::vector<double*> unknowns);
};

template<typename FlowFunction_t, typename AccumulationFunction_t>
    void NewtonRaphson<FlowFunction_t, AccumulationFunction_t>::solve(){
    _solver.compute(_jacobian);
    _solution_delta = _solver.solve(_residual);
};

template<typename FlowFunction_t, typename AccumulationFunction_t>
    void NewtonRaphson<FlowFunction_t, AccumulationFunction_t>::updateInteration(std::vector<double*> unknowns){
    
};
template<typename FlowFunction_t, typename AccumulationFunction_t>
    void NewtonRaphson<FlowFunction_t,AccumulationFunction_t>::iterate(const double relative_change_in_residual, const int term, const Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, Rock& rock){
    
    do{
        //Residual calculation
        for(auto fluid : characterized_fluids){
            for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){
                //calculateProperties(term, fluid->index(), cell->index());
                calculateResidual(term,cell);
            }
        }
        
    }while(_residual.squaredNorm()/_initial_residual.squaredNorm() > relative_change_in_residual);
    
};

#endif /* NEWTONRAPHSON_H */
