#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <Eigen/Sparse>
#include <cmath>

#include "Global.h"

#include "Mesh.h"
#include "Rock.h"
#include "Well.h"
#include "Equilibrium_Relation.h"

template<typename PropertiesFunction_t, typename FlowFunction_t, typename AccumulationFunction_t, typename PerforationFunction_t, typename WellFunction_t>
    class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat_t;
    typedef Eigen::Triplet<double> Tripletd_t;
    
    Eigen::BiCGSTAB<SparseMat_t, Eigen::IncompleteLUT<double, int> > _solver;

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
    PerforationFunction_t *_calculatePerforation;
    WellFunction_t *_calculateWellFlow;
    
    const inline int locate(const auto type, int input_selector, int input_index){
        using namespace Global;
        if(type == "fluid"){
            return cells_number*input_selector + input_index;
        }else{
            return cells_number*fluids_quantity + input_index;
        };
    };
    
 public:
    
 NewtonRaphson(PropertiesFunction_t calculateProperties, FlowFunction_t calculateFlow, AccumulationFunction_t calculateAccumulation, PerforationFunction_t calculatePerforation, WellFunction_t calculateWellFlow) : _calculateProperties(calculateProperties),
        _calculateFlow(calculateFlow), _calculateAccumulation(calculateAccumulation), _calculatePerforation(calculatePerforation), _calculateWellFlow(calculateWellFlow){};

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

    double inline calculateWellResidual(const int& term, std::shared_ptr<Well>& well){
        return well->operativeCondition()->value() - well->flow(term);
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
                
                residual_selector = equation->index();
                
                if(equation->type() == typeid(Well).name()){
                    
                    constexpr auto residual_type = "well";
                    auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(equation);
                    
                    row = locate(residual_type, residual_selector, residual_well->index());
                    _calculateWellFlow(term, residual_well);
                    _residual(row) = calculateWellResidual(term, residual_well);;
                    
                }else{
                    
                    constexpr auto residual_type = "fluid";
                    auto residual_fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(equation);
                    
                    residual_selector = residual_fluid->index();
            
                    for(auto cell = mesh.begin(); cell != mesh.end(); ++cell){
                
                        cell_index = cell->index();
                        _calculateProperties(term, *cell, rock);

                        row = locate(residual_type, residual_selector, cell_index);
                        
                        _residual(row) = calculateResidual(term,*residual_fluid, mesh, *cell, rock);
                        
                    };
                };
            };

            // This should be an equation component
            for(auto residual : equations){

                residual_selector = residual->index();
                
                if(residual->status()){
                
                    if(residual->type() == typeid(Well).name()){
                        
                        constexpr auto residual_type = "well";
                        auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(residual);
                        
                        for(auto variable : equations){

                            variable_selector = variable->index();
                            
                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){
                                    
                                    if(residual->index() == variable->index()){
                                        
                                        constexpr auto variable_type = "well";
                                        auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);

                                        row = locate(residual_type, residual_selector, residual_well->index());
                                        col = locate(variable_type, variable_selector, well_variable->index());

                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)+_machine_epsilon);
                                        _calculateWellFlow(term, well_variable);
                                        modified_residual = calculateWellResidual(term, well_variable);
                                    
                                    
                                        double derivative = (modified_residual - _residual(row))/_machine_epsilon;
                                        if(derivative != 0){
                                            _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                        };
                                                        
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)-_machine_epsilon);
                                        _calculateFlow(term, well_variable);
                                        
                                    };                     
                                    
                                }else{
                                    constexpr auto variable_type = "fluid";
                                    
                                    auto fluid_variable = std::dynamic_pointer_cast<Fluid,Equation_Base>(variable);
                                    for(auto perforation = residual_well->begin(); perforation!=residual_well->end(); ++perforation){
                                        auto cell = mesh.cell((*perforation)->index());

                                        auto unmodified_flow = (*perforation)->totalFlow();
                                        
                                        row = locate(residual_type, residual_selector, residual_well->index());
                                        col = locate(variable_type, variable_selector, cell.index());

                                        modified_epsilon = _machine_epsilon;

                                        modifyVariable(term, *fluid_variable, cell, modified_epsilon);
                                        
                                        _calculateProperties(term, cell, rock);

                                        _calculatePerforation(term, residual_well, *perforation);

                                        auto modified_flow = (*perforation)->totalFlow();

                                        double derivative = (modified_flow - unmodified_flow)/_machine_epsilon;
                                        if(derivative != 0){
                                            _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                        };
                                        
                                        modified_epsilon = -_machine_epsilon;

                                        modifyVariable(term, *fluid_variable, cell, modified_epsilon);
                                        
                                        _calculateProperties(term, cell, rock);

                                        _calculatePerforation(term, residual_well, *perforation);
                                        
                                    };
                                    
                                };
                                
                            };
                            
                        };
                        
                    }else{
                        
                        constexpr auto residual_type = "fluid";
                        auto residual_fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(residual);
                
                        residual_selector = residual_fluid->index();

                        //This should be a principal variable 
                        for(auto variable : equations){
                            
                            variable_selector = variable->index();
                            
                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){
                                    
                                    constexpr auto variable_type = "well";
                                    
                                    auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);
                                    for(auto perforation = well_variable->begin(); perforation!=well_variable->end(); ++perforation){

                                        auto cell = mesh.cell((*perforation)->index());
                                        
                                        row = locate(residual_type, residual_selector, cell.index());
                                        col = locate(variable_type, variable_selector, well_variable->index());
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)+_machine_epsilon);
                                        _calculatePerforation(term, well_variable, *perforation);
                                    
                                        modified_residual = calculateResidual(term,*residual_fluid, mesh,cell,rock);
                                        
                                        double derivative = (modified_residual - _residual(row))/_machine_epsilon;
                                        if(derivative != 0){
                                            _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                        };
                                        
                                    
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)-_machine_epsilon);
                                        _calculatePerforation(term, well_variable, *perforation);
                                    };
                                    
                                }else{

                                    constexpr auto variable_type = "fluid";
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
                                            row = locate(residual_type, residual_selector, neighbor_index);
                                            col = locate(variable_type, variable_selector, cell_index);
                                            modified_residual = calculateResidual(term,*residual_fluid, mesh,*neighbor_cell,rock);

                                            double derivative = (modified_residual - _residual(neighbor_index))/_machine_epsilon;
                                            if(derivative != 0){
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };
                        
                                        int row = locate(residual_type, residual_selector, cell_index);
                                        int col = locate(variable_type, variable_selector, cell_index);
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
            
            if(equation->status()){

                if(equation->type() == typeid(Well).name()){
                    
                    constexpr auto residual_type = "well";
                    auto well = std::dynamic_pointer_cast<Well,Equation_Base>(equation);
                    row = locate(residual_type, residual_selector, well->index());

                    well->boreholePressure(term, well->boreholePressure(term)+_solution_delta(row));
                    
                }else{

                    constexpr auto residual_type = "fluid";
                    auto fluid = std::dynamic_pointer_cast<Fluid,Equation_Base>(equation);
            
                    residual_selector = fluid->index();
            
                    for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){

                        cell_index = cell->index();
                        row = locate(residual_type,residual_selector,cell_index);

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
