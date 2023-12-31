#ifndef NEWTONRAPHSON_H
#define NEWTONRAPHSON_H

#include <Eigen/Sparse>
#include <cmath>

#include "Initial_Conditions.h"
#include "Database.h"

#include "Mesh.h"
#include "Rock.h"
#include "Production_Well.h"
#include "Injection_Well.h"

template<typename PropertiesFunction_t, typename FlowFunction_t, typename AccumulationFunction_t, typename PerforationFunction_t, typename WellFunction_t, typename EstimatorFunction_t>
    class NewtonRaphson{
 private:
    
    typedef Eigen::SparseMatrix<double> SparseMat_t;
    typedef Eigen::Triplet<double> Tripletd_t;
    
    Eigen::BiCGSTAB<SparseMat_t, Eigen::IncompleteLUT<double, int> > _solver;

    const double _machine_epsilon = 2.0*std::sqrt(std::numeric_limits<double>::epsilon());
    const double _relative_change_in_residual=1e-11;
    const double _divergence_tolerance = 1e4;

    std::vector<Tripletd_t> _non_zeros;
    
    SparseMat_t _jacobian;
    Eigen::VectorXd _residual;
    Eigen::VectorXd _past_residual;
    Eigen::VectorXd _initial_residual;
    Eigen::VectorXd _solution_delta;

    double _aux_variable;

    PropertiesFunction_t *_calculateProperties;
    FlowFunction_t *_calculateFlow;
    AccumulationFunction_t *_calculateAccumulation;
    PerforationFunction_t *_calculatePerforation;
    WellFunction_t *_calculateWellFlow;
    EstimatorFunction_t* _estimatePressure;
    
    const inline int locate(const auto type, int input_selector, int input_index){
        using namespace Initial_Conditions;
        if(type == "phase"){
            return cells_quantity*input_selector + input_index;
        }else{
            return cells_quantity*phases_quantity + input_index;
        };
    };
    
 public:

    NewtonRaphson(const NewtonRaphson& other_newton){
        _jacobian = other_newton._jacobian;
        _residual = other_newton._residual;
        _initial_residual = other_newton._initial_residual;
        _solution_delta = other_newton._solution_delta;

        _non_zeros = other_newton._non_zeros;

        _calculateProperties = other_newton._calculateProperties;
        _calculateFlow = other_newton._calculateFlow;
        _calculateAccumulation = other_newton._calculateAccumulation;
        _calculatePerforation = other_newton._calculatePerforation;
        _calculateWellFlow = other_newton._calculateWellFlow;
        _estimatePressure = other_newton._estimatePressure;
    };
    
 NewtonRaphson(const int& mat_size,const int& max_number_of_non_zeros, PropertiesFunction_t calculateProperties, FlowFunction_t calculateFlow, AccumulationFunction_t calculateAccumulation, PerforationFunction_t calculatePerforation, WellFunction_t calculateWellFlow, EstimatorFunction_t estimatePressure) : _calculateProperties(calculateProperties),
        _calculateFlow(calculateFlow), _calculateAccumulation(calculateAccumulation), _calculatePerforation(calculatePerforation), _calculateWellFlow(calculateWellFlow), _estimatePressure(estimatePressure){

        //= phases_quantity*cells_quantity + wells_with_equation;

        // = phases_quantity*phases_quantity*cells_quantity + total_of_perforations*wells_with_equation;
        
        _jacobian.resize(mat_size, mat_size);

        _jacobian.reserve (max_number_of_non_zeros);
        _non_zeros.reserve(max_number_of_non_zeros);

        _residual         = Eigen::VectorXd(mat_size);
        _initial_residual = Eigen::VectorXd(mat_size);
        _solution_delta   = Eigen::VectorXd(mat_size);

        _residual.setZero();
        _initial_residual.setZero();
        _solution_delta.setZero();
    };

    void modifyVariable(const int& term, Phase& phase, const std::shared_ptr<Cell>& cell, double& modified_epsilon, const bool undo){


        double scaled_epsilon = 0.0;
        
        const int cell_index = cell->index();
            
        if(phase.main()){
            if(!undo){
                _aux_variable = phase.pressure(term, cell_index);
                
                scaled_epsilon = std::abs(_aux_variable) * _machine_epsilon;
                if(scaled_epsilon == 0.0) scaled_epsilon = _machine_epsilon;
                
                phase.pressure(term, cell_index, phase.pressure(term, cell_index)+scaled_epsilon);
                
                modified_epsilon = phase.pressure(term, cell_index)-_aux_variable;
                //std::cout << modified_epsilon;
            }else{
                phase.pressure(term, cell_index, _aux_variable);
            };
        }else{
            if(!undo){
                _aux_variable = phase.saturation(term, cell_index);
                
                scaled_epsilon = std::abs(_aux_variable) * _machine_epsilon;
                if(scaled_epsilon == 0.0) scaled_epsilon = _machine_epsilon;
                
                phase.saturation(term, cell_index, phase.saturation(term, cell_index)+scaled_epsilon);
                
                modified_epsilon = phase.saturation(term, cell_index)-_aux_variable;
                
            }else{
                phase.saturation(term, cell_index, _aux_variable);
            };
        };
    };
    
    double calculateResidual(const int& term, Phase& phase, const Mesh& mesh, const std::shared_ptr<Cell>& cell, Rock& rock, std::vector<std::shared_ptr<Well>>& wells)
    {
        double flow=0.0;
        for (auto face = cell->begin(); face!=cell->end(); ++face){
            flow = flow + _calculateFlow(term, phase, mesh, cell, *face, rock);
        };
        double accumulation = _calculateAccumulation(term, phase, cell, rock);

        double well_contribution = 0.0;
        
        for (auto well : wells){
            for(auto perforation = well->begin(); perforation !=well->end(); ++perforation){
                if((*perforation)->index() == cell->index()){
                    if((*perforation)->type() == typeid(Production_Perforation).name()){
                        auto production_perf = std::dynamic_pointer_cast<Production_Perforation, Perforation>(*perforation);
                        well_contribution += production_perf->flow(phase.index());
                        //std::cout <<"\n"<< well_contribution;
                    }else{
                        auto injection_well = std::dynamic_pointer_cast<Injection_Well, Well>(well);
                        if (injection_well->injectionPhase()->index() == phase.index()){
                            well_contribution += (*perforation)->totalFlow();
                            
                        };
                    };
                };
            };
        };

        return accumulation - flow - well_contribution;
    };

    double inline calculateWellResidual(const int& term, std::shared_ptr<Well>& well){
        //std::cout << well->operativeCondition()->value() - well->flow(term);
        return well->operativeCondition()->value() - well->flow(term);
    };

    void solve(){

        _residual = -_residual;
         //std::cout << "Jacobian: \n"<< _jacobian<<std::endl;
        _solver.preconditioner().setDroptol(0.00000001);
        _solver.compute(_jacobian);
        _solution_delta = _solver.solve(_residual);
        //std::cout << "Solution Delta: "<< _solution_delta<<std::endl;
        _residual = -_residual;
    };
    
    double iterate(const int& term, const Mesh& mesh, std::vector<std::shared_ptr<Well>>& wells, std::vector<std::shared_ptr<Equation_Base>>& equations, Rock& rock){

        int residual_selector;
        int cell_index;
        int variable_selector;
        double modified_epsilon;
        double scaled_epsilon;
        double modified_residual;
        int row;
        int col;
        int iteration=0;
        double stability_tolerance=1e4;
        double initial_norm;
        //decltype(_residual) past_residual = _residual;

        double tolerance=0;
        double past_tolerance=0;
        
    
        do{

            _residual.setZero();
            _solution_delta.setZero();

            for(auto cell = mesh.begin(); cell != mesh.end(); ++cell){
                _calculateProperties(term, *cell, rock);
            };

            for(auto well : wells){
                
                if(well->operativeStatus()==1){

                    if(well->operativeCondition()->type() == "FLOW"){
                        _estimatePressure(term, well);
                    };
                    
                    well->operativeStatus(0);
                };
                _calculateWellFlow(term, well);
            };
            
            //Residual calculation
            residual_selector = 0;
            
            for(auto equation : equations){

                if(equation->status()){
                
                    if(equation->type() == typeid(Well).name()){
                    
                        constexpr auto residual_type = "well";
                        auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(equation);
                    
                        row = locate(residual_type, residual_selector, residual_well->index());
                        _residual(row) = calculateWellResidual(term, residual_well);
                    
                    }else{
                    
                        constexpr auto residual_type = "phase";
                        auto residual_phase = std::dynamic_pointer_cast<Phase,Equation_Base>(equation);
            
                        for(auto cell = mesh.begin(); cell != mesh.end(); ++cell){
                
                            cell_index = (*cell)->index();

                            row = locate(residual_type, residual_selector, cell_index);
                        
                            _residual(row) = calculateResidual(term,*residual_phase, mesh, *cell, rock, wells);
                        
                        };
                    };
                
                    ++residual_selector;
                    
                };
            };

            //std::cout << "Residual: "<< _residual << std::endl;
            
            //Jacobian Calculation

            residual_selector = 0;
            // This should be an equation component
            for(auto residual : equations){

                //residual_selector = residual->index();
                variable_selector = 0;
                if(residual->status()){
                
                    if(residual->type() == typeid(Well).name()){
                        
                        constexpr auto residual_type = "well";
                        auto residual_well = std::dynamic_pointer_cast<Well,Equation_Base>(residual);

                        
                        for(auto variable : equations){

                            //variable_selector = variable->index();
                            
                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){
                                    
                                    if(residual->index() == variable->index()){
                                        
                                        constexpr auto variable_type = "well";
                                        auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);
                                        scaled_epsilon=0.0;

                                        row = locate(residual_type, residual_selector, residual_well->index());
                                        col = locate(variable_type, variable_selector, well_variable->index());
                                        
                                        _aux_variable = well_variable->boreholePressure(term);
                                        scaled_epsilon = std::abs(_aux_variable) * _machine_epsilon;
                                        
                                        if(scaled_epsilon == 0.0) {scaled_epsilon = _machine_epsilon;};
                                        
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)+scaled_epsilon);

                                        modified_epsilon = well_variable->boreholePressure(term) - _aux_variable;

                                        //std::cout << modified_epsilon << std::endl;
                                            
                                        _calculateWellFlow(term, well_variable);
                                        modified_residual = calculateWellResidual(term, well_variable);
                                    
                                    
                                        double derivative = (modified_residual - _residual(row))/modified_epsilon;

                                        if(std::isnan(derivative)){derivative = 0.0;}; 
                                        
                                        if(derivative != 0.0 || row == col){
                                            if(derivative == 0 && row == col){
                                                _non_zeros.push_back(Tripletd_t(row,col,_machine_epsilon));
                                            }else{
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };
                                                        
                                        well_variable->boreholePressure(term, _aux_variable);
                                        _calculateWellFlow(term, well_variable);
                                        
                                    };                     
                                    
                                }else{
                                    constexpr auto variable_type = "phase";
                                    
                                    auto phase_variable = std::dynamic_pointer_cast<Phase,Equation_Base>(variable);
                                    for(auto perforation = residual_well->begin(); perforation!=residual_well->end(); ++perforation){
                                        auto cell = mesh.cell((*perforation)->index());

                                        scaled_epsilon=0.0;

                                        auto unmodified_flow = (*perforation)->totalFlow();
                                        
                                        row = locate(residual_type, residual_selector, residual_well->index());
                                        col = locate(variable_type, variable_selector, cell->index());

                                        modified_epsilon = _machine_epsilon;

                                        modifyVariable(term, *phase_variable, cell, modified_epsilon, false);
                                        
                                        _calculateProperties(term, cell, rock);

                                        _calculatePerforation(term, residual_well, *perforation);

                                        auto modified_flow = (*perforation)->totalFlow();

                                        //std::cout << modified_epsilon << std::endl;

                                        double derivative = (modified_flow - unmodified_flow)/modified_epsilon;
                                        if(std::isnan(derivative)){derivative = 0.0;};
                                        
                                        if(derivative != 0.0 || row == col){
                                            if(derivative == 0 && row == col){
                                                _non_zeros.push_back(Tripletd_t(row,col,_machine_epsilon));
                                            }else{
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };
                                        
                                        modified_epsilon = -_machine_epsilon;

                                        modifyVariable(term, *phase_variable, cell, modified_epsilon, true);
                                        
                                        _calculateProperties(term, cell, rock);

                                        _calculatePerforation(term, residual_well, *perforation);
                                        
                                    };
                                    
                                };
                                
                                ++variable_selector;
                                
                            };
                            
                        };
                        
                    }else{
                        
                        constexpr auto residual_type = "phase";
                        auto residual_phase = std::dynamic_pointer_cast<Phase,Equation_Base>(residual);
                
                        //residual_selector = residual_phase->index();

                        //This should be a main variable 
                        for(auto variable : equations){
                            
                            //variable_selector = variable->index();
                            
                            if(variable->status()){

                                if(variable->type() == typeid(Well).name()){
                                    
                                    constexpr auto variable_type = "well";
                                    
                                    auto well_variable = std::dynamic_pointer_cast<Well,Equation_Base>(variable);
                                    for(auto perforation = well_variable->begin(); perforation!=well_variable->end(); ++perforation){

                                        auto cell = mesh.cell((*perforation)->index());

                                        scaled_epsilon=0.0;
                                        
                                        row = locate(residual_type, residual_selector, cell->index());
                                        col = locate(variable_type, variable_selector, well_variable->index());
                                        
                                        _aux_variable = well_variable->boreholePressure(term);
                                        
                                        scaled_epsilon = std::abs(_aux_variable) * _machine_epsilon;
                                        if(scaled_epsilon == 0.0) scaled_epsilon = _machine_epsilon;
                                        
                                        well_variable->boreholePressure(term, well_variable->boreholePressure(term)+scaled_epsilon);
                                        
                                        //std::cout << well_variable->boreholePressure(term);
                                        
                                        modified_epsilon = well_variable->boreholePressure(term) - _aux_variable;
                                        
                                        _calculatePerforation(term, well_variable, *perforation);
                                    
                                        modified_residual = calculateResidual(term,*residual_phase, mesh,cell,rock, wells);

                                        //std::cout << modified_epsilon << std::endl;
                                        
                                        double derivative = (modified_residual - _residual(row))/modified_epsilon;
                                        if(std::isnan(derivative)){derivative = 0.0;};
                                        
                                        if(derivative != 0.0 || row == col){
                                            if(derivative == 0 && row == col){
                                                _non_zeros.push_back(Tripletd_t(row,col,_machine_epsilon));
                                            }else{
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };
                                        
                                    
                                        well_variable->boreholePressure(term, _aux_variable);
                                        _calculatePerforation(term, well_variable, *perforation);
                                    };
                                    
                                }else{

                                    constexpr auto variable_type = "phase";
                                    auto phase_variable = std::dynamic_pointer_cast<Phase,Equation_Base>(variable);
                                    variable_selector = phase_variable->index();
                    
                                    for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){

                                        modified_epsilon = _machine_epsilon;
                                        cell_index = (*cell)->index();
                        
                                        modifyVariable(term, *phase_variable, *cell, modified_epsilon, false);
                                        _calculateProperties(term, *cell, rock);
                        
                                        for (auto face = (*cell)->begin(); face!=(*cell)->end(); ++face){
                            
                                            auto neighbor_cell = (*face)->neighborCell().lock();
                                            int neighbor_index = neighbor_cell->index();
                                            row = locate(residual_type, residual_selector, neighbor_index);
                                            col = locate(variable_type, variable_selector, cell_index);
                                            modified_residual = calculateResidual(term,*residual_phase, mesh,neighbor_cell,rock, wells);

                                            //std::cout << modified_epsilon << std::endl;
                                            
                                            double derivative = (modified_residual - _residual(neighbor_index))/modified_epsilon;
                                            if(std::isnan(derivative)){derivative = 0.0;};
                                            
                                            if(derivative != 0 || row == col){
                                                if(derivative == 0 && row == col){
                                                    _non_zeros.push_back(Tripletd_t(row,col,_machine_epsilon));
                                                }else{
                                                    _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                                };
                                            };

                                            neighbor_cell.reset();
                                        };
                        
                                        int row = locate(residual_type, residual_selector, cell_index);
                                        int col = locate(variable_type, variable_selector, cell_index);
                                        modified_residual = calculateResidual(term,*residual_phase, mesh, *cell, rock, wells);

                                        //std::cout << modified_epsilon << std::endl;

                                        double derivative = (modified_residual - _residual(cell_index))/modified_epsilon;
                                        if(std::isnan(derivative)){derivative = 0.0;};
                                        
                                        if(derivative != 0.0 || row == col){
                                            if(derivative == 0 && row == col){
                                                _non_zeros.push_back(Tripletd_t(row,col,_machine_epsilon));
                                            }else{
                                                _non_zeros.push_back(Tripletd_t(row,col,derivative));
                                            };
                                        };

                                        modified_epsilon = -_machine_epsilon;
                                        modifyVariable(term, *phase_variable, *cell, modified_epsilon, true);
                                        _calculateProperties(term, *cell, rock);
                                    
                                    };
                                    
                                };
                                
                                ++variable_selector;
                              
                            };
                            
                        };
                        
                    };

                    ++residual_selector;
                };
                
            };

            _jacobian.setFromTriplets(_non_zeros.begin(), _non_zeros.end());

            if(iteration == 0){
                _initial_residual = _residual;
            }
            
            solve();
            update(term, mesh, equations);
            _jacobian.setZero();
            _non_zeros.clear();

            //auto difference = past_residual - _residual;
            //
            //tolerance = difference.squaredNorm()/_residual.squaredNorm();

            initial_norm = _initial_residual.squaredNorm();
            
            tolerance = _residual.squaredNorm();

            if(iteration > 0){

                past_tolerance = _past_residual.squaredNorm();
                
                stability_tolerance = std::abs(past_tolerance/tolerance) - 1.0;
            };
 
            _past_residual = _residual;

            //std::cout << "Tolerance: " << tolerance << " at iteration: " << iteration << std::endl;
            ++iteration;
            
        }while(tolerance > _relative_change_in_residual && tolerance < _divergence_tolerance && iteration < 100);

        return tolerance;
    };
    
    void update(const int& term, const Mesh& mesh, std::vector<std::shared_ptr<Equation_Base>>& equations){
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

                    //Well at atmospheric pressure [pa]
                    if(well->boreholePressure(term) < 101325.0){
                        well->boreholePressure(term, 101325.0);
                        well->operativeCondition()->type("SHUT");
                    }
                    
                }else{

                    constexpr auto residual_type = "phase";
                    auto phase = std::dynamic_pointer_cast<Phase,Equation_Base>(equation);
            
                    residual_selector = phase->index();
            
                    for(auto cell = mesh.begin(); cell !=mesh.end(); ++cell){

                        cell_index = (*cell)->index();
                        row = locate(residual_type,residual_selector,cell_index);

                        if(phase->main()){
                            phase->pressure(term, cell_index, phase->pressure(term, cell_index)+_solution_delta(row));
                            if(phase->pressure(term, cell_index) < 101325.0){
                                phase->pressure(term, cell_index, 101325.0);
                            };
                            
                        }else{
                            phase->saturation(term, cell_index, phase->saturation(term, cell_index)+_solution_delta(row));
                            if(phase->saturation(term, cell_index) < 0){
                                phase->saturation(term, cell_index, 0);
                            };
                        };
                    };
                };
            };
        };
    };
};

#endif /* NEWTONRAPHSON_H */
