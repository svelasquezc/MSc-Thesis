#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <Eigen/Sparse>
#include <cmath>

#include "Mesh.h"
#include "Rock.h"
#include "Well.h"
#include "Equilibrium_Relation.h"

template<typename PropertiesFunction_t, typename FlowFunction_t, typename AccumulationFunction_t>
    class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat_t;
    typedef Eigen::Triplet<double> Tripletd_t;
    
    Eigen::BiCGSTAB<SparseMat_t, Eigen::IncompleteLUT<double, int> > _solver;

    mutable int _cells_number=0;

    int _iteration=0;

    const double _machine_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
    const double _relative_change_in_residual=1e-4;

    std::vector<Tripletd_t> _non_zeros;
    
    SparseMat_t _jacobian;
    Eigen::VectorXd _residual;
    Eigen::VectorXd _initial_residual;
    Eigen::VectorXd _solution_delta;

    PropertiesFunction_t *_calculateProperties;
    FlowFunction_t *_calculateFlow;
    AccumulationFunction_t *_calculateAccumulation;
    
    const inline int locate(int input_selector, int input_index){
        return _cells_number*input_selector + input_index;
    };
    
 public:
    
 NewtonRaphson(PropertiesFunction_t calculateProperties, FlowFunction_t calculateFlow, AccumulationFunction_t calculateAccumulation) : _calculateProperties(calculateProperties),
        _calculateFlow(calculateFlow), _calculateAccumulation(calculateAccumulation){};

    void modifyVariable(const int& term, Fluid& fluid, Cell& cell, const double modified_epsilon){
        
        const int cell_index = cell.index();
            
        if(fluid.principal()){
            fluid.pressure(term, cell_index, fluid.pressure(term, cell_index)+modified_epsilon);
        }else{
            fluid.saturation(term, cell_index, fluid.saturation(term, cell_index)+modified_epsilon);
        };
    };
    
    double calculateResidual(const int& term, Fluid& fluid, Mesh& mesh, Cell& cell, Rock& rock)
    {
        double flow=0;
        for (auto face = cell.begin(); face!=cell.end(); ++face){
            flow = flow + _calculateFlow(term, fluid, mesh, cell, *face, rock);
        };
        double accumulation = _calculateAccumulation(term, fluid, cell, rock);

        return accumulation - flow;
    };

    void solve(){
        _solver.compute(_jacobian);
        _solution_delta = _solver.solve(_residual);
    }
    
    void iterate(const int term, const Mesh& mesh, std::vector<std::shared_ptr<Equation_Base>>& equations, Rock& rock){

        int residual_selector;
        int cell_index;
        int variable_selector;
        double modified_epsilon;
        double modified_residual;
        int row;
        int col;
    
        do{
            //Residual calculation
            for(auto equation : equations){

                if(equation->type() == typeid(Well).name()){
                    auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(equation);
                }else{

                    auto residual_fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(equation);
                    
                    residual_selector = residual_fluid->index();
            
                    for(auto cell = mesh.begin(); cell != mesh.end(); ++cell){
                
                        cell_index = cell->index();
                        _calculateProperties(term, *cell, rock);
                        residual(locate(residual_selector, cell_index)) = calculateResidual(term,*residual_fluid, mesh, *cell, rock);
                        
                    };
                };
            };

            // This should be an equation component
            for(auto residual : equations){

                if(residual->status()){
                
                    if(residual->type() == typeid(Well).name()){

                        auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(residual);
                        
                        for(auto variable : equations){

                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){
                                    //int row = locate(residual_selector, well_index);
                                    //int col = locate(variable_selector, cell_index);
                                    
                                    auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);
                                    well_variable->boreholePressure(term, well_variable->boreholePressure(term)+_machine_epsilon);
                                    //calculateFlow(term, well_variable);
                                    //modified_residual = calculateWellResidual(well_variable);
                                    
                                    /*
                                      double derivative = (modified_residual - _residual(neighbor_index))/_machine_epsilon;
                                      if(derivative != 0){
                                      _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                      };
                                    */
                                    
                                    
                                    well_variable->boreholePressure(term, well_variable->boreholePressure(term)-_machine_epsilon);
                                    //calculateFlow(term, well_variable);
                     
                                    
                                }else{
                                    
                                    auto fluid_variable = std::dynamic_pointer_cast<Fluid,Equation_Base>(variable);
                                    for(auto perforation = residual_well->begin(); perforation!=residual_well->end(); ++perforation){
                                        auto cell = mesh.cell((*perforation)->index());

                                        auto unmodified_flow = (*perforation)->totalFlow();
                                        
                                        //int col = locate(variable_selector, cellindex());

                                        modified_epsilon = _machine_epsilon;

                                        modifyVariable(term, *fluid_variable, cell, modified_epsilon);
                                        
                                        _calculateProperties(term, cell, rock);

                                        //calculatePerforation(term, well_variable, *perforation);

                                        auto modified_flow = (*perforation)->totalFlow();

                                        double derivative = (modified_flow - unmodified_flow)/_machine_epsilon;
                                        if(derivative != 0){
                                            _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                        };
                                        
                                        modified_epsilon = -_machine_epsilon;

                                        modifyVariable(term, *fluid_variable, cell, modified_epsilon);
                                        
                                        _calculateProperties(term, cell, rock);

                                        //calculatePerforation(term, well_variable, *perforation);
                                        
                                    };
                                    
                                };
                                
                            };
                            
                        };
                        
                    }else{

                        auto residual_fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(residual);
                
                        residual_selector = residual_fluid->index();

                        //This should be a principal variable 
                        for(auto variable : equations){

                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){

                                    auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);
                                    for(auto perforation = well_variable->begin(); perforation!=well_variable->end(); ++perforation){

                                        auto cell = mesh.cell((*perforation)->index());
                                        
                                        //int row = locate(residual_selector, neighbor_index);
                                        //int col = locate(variable_selector, cellindex());
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)+_machine_epsilon);
                                        //calculatePerforation(term, well_variable, *perforation);
                                    
                                        modified_residual = calculateResidual(term,*residual_fluid, mesh,cell,rock);
                                    
                                        /*
                                          double derivative = (modified_residual - _residual(row))/_machine_epsilon;
                                          if(derivative != 0){
                                          _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                          };
                                        */
                                    
                                    
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)-_machine_epsilon);
                                        //calculatePerforation(term, well_variable, *perforation);
                                    };
                                }else{

                                    auto fluid_variable = std::dynamic_pointer_cast<Fluid,Equation_Base>(variable);
                                    variable_selector = fluid_variable->index();
                    
                                    for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){

                                        modified_epsilon = _machine_epsilon;
                                        cell_index = cell->index();
                        
                                        modifyVariable(term, *fluid_variable, *cell, modified_epsilon);
                                        _calculateProperties(term, *cell, rock);
                        
                                        for (auto face = cell->begin(); face!=cell->end(); ++face){
                            
                                            auto neighbor_cell = face->neighbor();
                                            int neighbor_index = neighbor_cell->index();
                                            row = locate(residual_selector, neighbor_index);
                                            col = locate(variable_selector, cell_index);
                                            modified_residual = calculateResidual(term,*residual_fluid, mesh,*neighbor_cell,rock);

                                            double derivative = (modified_residual - _residual(neighbor_index))/_machine_epsilon;
                                            if(derivative != 0){
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };
                        
                                        int row = locate(residual_selector, cell_index);
                                        int col = locate(variable_selector, cell_index);
                                        modified_residual = calculateResidual(term,*residual_fluid, mesh, *cell, rock);

                                        double derivative = (modified_residual - _residual(cell_index))/_machine_epsilon;
                                        if(derivative != 0){
                                            _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                        };

                                        modified_epsilon = -_machine_epsilon;
                                        modifyVariable(term, *fluid_variable, *cell, modified_epsilon);
                                        _calculateProperties(term, *cell, rock);
                                    
                                    };
                                    
                                };
                                
                            };
                            
                        };
                        
                    };
                    
                };
                
            };

            _jacobian.setFromTriplets(_non_zeros.begin(), _non_zeros.end());
            solve();
            update(term, mesh, equations);
            _non_zeros.clear();
        
        }while(_residual.squaredNorm()/_initial_residual.squaredNorm() > _relative_change_in_residual);
    
    };
    
    void update(const int term, const Mesh& mesh, std::vector<std::shared_ptr<Equation_Base>>& equations){
        int cell_index;
        int residual_selector;
        int row;

        for(auto equation : equations){
            
            if(variable->status()){

                if(variable->type() == typeid(Well).name()){

                    auto well = std::dynamic_pointer_cast<Well,Equation_Base>(equation);
                    
                }else{
                    
                    auto fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(equation);
            
                    residual_selector = fluid->index();
            
                    for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){

                        cell_index = cell->index();
                        row = locate(residual_selector,cell_index);

                        if(fluid->principal()){
                            fluid->pressure(term, cell_index, fluid->pressure(term, cell_index)+_solution_delta(row));
                        }else{
                            fluid->saturation(term, cell_index, fluid->saturation(term, cell_index)+_solution_delta(row));
                        };
                    };
                };
            };
        };
    };
};

#endif /* NEWTONRAPHSON_H */
