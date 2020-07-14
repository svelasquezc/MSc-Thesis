#ifndef WELL_FUNCTIONS_H
#define WELL_FUNCTIONS_H

#include "Production_Well.h"
#include "Injection_Well.h"

void calculateGeometry(const int& term, std::shared_ptr<Well>& well, std::shared_ptr<Perforation>& perforation){
    
    auto cell = defined_mesh->cell(perforation->index());

    auto numeration = cell->numeration3D();

    auto absolute_permeability= characterized_rock->absolutePermeability(term,cell->index());
    
    int axis_x = numeration[0];
    int axis_y = numeration[1];
    int axis_z = numeration[2];

    perforation->calculateEquivalentRadius(term, absolute_permeability[1],
                                           absolute_permeability[0],
                                           defined_mesh->thickness(1,axis_y),
                                           defined_mesh->thickness(0,axis_x)
                                           );
    
    perforation->calculateWellIndex       (term,
                                           defined_mesh->thickness(2,axis_z),
                                           absolute_permeability[1],
                                           absolute_permeability[0],
                                           well->radius()
                                           );

    
};

void estimateWellPressure(const int& term, std::shared_ptr<Well>& well){

    double perforation_mobility = 0;
    double hydrostatic_head = 0;
    double accumulated_mobility = 0;
    double estimated_flow=0;
    double accumulated_flow=0;
    
    
    std::shared_ptr<Injection_Perforation> injection_perf;
    std::shared_ptr<Production_Perforation> production_perf;

    for(auto perforation = well->begin(); perforation !=well->end(); ++perforation){    
        calculateGeometry(term, well, *perforation);
    };

    std::shared_ptr<Injection_Well> injection_well;
    std::shared_ptr<Production_Well> production_well;

    std::shared_ptr<Phase> main_phase;

    for(auto phase : characterized_phases){
        if(phase->main()){
            main_phase = phase;
        };
    }
    
    if(well->type() == typeid(Production_Well).name()){
        production_well = std::dynamic_pointer_cast<Production_Well, Well>(well);
        for(auto perforation = production_well->begin(); perforation !=production_well->end(); ++perforation){
            production_perf = std::dynamic_pointer_cast<Production_Perforation, Perforation>(*perforation);

            auto cell = defined_mesh->cell(production_perf->index());
            
            for(auto phase : characterized_phases){
                
                if(phase->type() != "Gas"){
                    perforation_mobility = production_perf->wellIndex(term)*
                        phase->relativePermeability(term, cell->index())/(phase->formationVolumeFactor(term, cell->index())*phase->viscosity(term, cell->index()));

                    accumulated_mobility += perforation_mobility;

                    hydrostatic_head = phase->density(term, cell->index())*gravity*(production_well->boreholeDepth()-cell->depth());

                    estimated_flow = perforation_mobility*(main_phase->pressure(term, cell->index()) + hydrostatic_head);
                    accumulated_flow += estimated_flow;
                    
                };
                
            };
            
        };
        
    }else{
        
        injection_well = std::dynamic_pointer_cast<Injection_Well, Well>(well);
        auto injection_phase = injection_well->injectionPhase();
        
        for(auto perforation = injection_well->begin(); perforation !=injection_well->end(); ++perforation){
            injection_perf = std::dynamic_pointer_cast<Injection_Perforation, Perforation>(*perforation);

            auto cell = defined_mesh->cell(injection_perf->index());

            double total_mobility = 0;

            std::shared_ptr<Phase> main_phase;
        
            for(auto phase : characterized_phases){
                if(phase->main()){
                    main_phase = phase;
                };
                total_mobility += phase->relativePermeability(term, cell->index())/phase->viscosity(term, cell->index());
            };
            
            perforation_mobility = injection_perf->wellIndex(term)*total_mobility
                /injection_phase->formationVolumeFactor(term, cell->index());

            accumulated_mobility += perforation_mobility;

            hydrostatic_head = injection_phase->density(term, cell->index())*gravity*(injection_well->boreholeDepth()-cell->depth());

            estimated_flow = perforation_mobility*(main_phase->pressure(term, cell->index()) + hydrostatic_head);
            accumulated_flow += estimated_flow;
                    
        };
        
    };

    well->boreholePressure(term, (accumulated_flow + well->flow(term))/accumulated_mobility);
};

double calculatePeacemanProduction(const int& term, const double main_pressure, const std::shared_ptr<Phase>& phase, const std::shared_ptr<Cell>& cell, const double well_index, const double borehole_pressure, const double borehole_depth){
    
    auto cell_index = cell->index();
    double peaceman_flow = 0.0;
    
    peaceman_flow = (well_index * phase->relativePermeability(term, cell_index) / ((phase->formationVolumeFactor(term, cell_index)*phase->viscosity(term, cell_index)))) *
        (borehole_pressure - main_pressure - (phase->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));
    
    //std::cout << "\n mira::wi  " <<    well_index ;
    //std::cout << "\n mira::kr  " << phase->relativePermeability(term, cell_index);
    //std::cout << "\n mira::bol  " << phase->formationVolumeFactor(term, cell_index);
    //std::cout << "\n mira::vis  " << phase->viscosity(term, cell_index);
    //std::cout << "\n mira::bhp  " << borehole_pressure;
    //std::cout << "\n mira::bhp  " <<    main_pressure;
    //std::cout << "\n mira::Pres  " << (borehole_pressure - main_pressure - (phase->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));

    
    //std::cout <<"\n mira:" << borehole_depth - cell->depth() << "\n";

    return peaceman_flow;
    
};

double calculatePeacemanInjection(const int& term, const double main_pressure, const double total_mobility, const std::shared_ptr<Phase>& phase, const std::shared_ptr<Cell>& cell, const double well_index, const double borehole_pressure, const double borehole_depth){
    
    auto cell_index = cell->index();
    double peaceman_flow = 0;

    peaceman_flow = (well_index * total_mobility / phase->formationVolumeFactor(term, cell_index)) *
        (borehole_pressure - main_pressure - (phase->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));


    return peaceman_flow;
};

void calculatePerforation(const int& term, std::shared_ptr<Well>& well, std::shared_ptr<Perforation>& perforation){

    std::shared_ptr<Injection_Perforation> injection_perf;
    std::shared_ptr<Production_Perforation> production_perf;

    
    std::shared_ptr<Injection_Well> injection_well;
    std::shared_ptr<Production_Well> production_well;

    if(well->type() == typeid(Production_Well).name()){
        production_well = std::dynamic_pointer_cast<Production_Well, Well>(well);
    }else{
        injection_well = std::dynamic_pointer_cast<Injection_Well, Well>(well);
    };

    double main_pressure;

    for(auto phase : characterized_phases){
        if(phase->main()){
            main_pressure = phase->pressure(term, perforation->index());
        }
    }

    calculateGeometry(term, well, perforation);
    
    if(perforation->type() == typeid(Production_Perforation).name()){
        
        production_perf = std::dynamic_pointer_cast<Production_Perforation, Perforation>(perforation);
        
        for(auto phase : characterized_phases){
            production_perf->flow(phase->index(),
                                  calculatePeacemanProduction(term,
                                                              main_pressure,
                                                              phase,
                                                              defined_mesh->cell(production_perf->index()),
                                                              production_perf->wellIndex(term),
                                                              well->boreholePressure(term),
                                                              well->boreholeDepth()
                                                              ));
        };
        
    }else{
        
        double total_mobility = 0;
        
        for(auto phase : characterized_phases){
            total_mobility += phase->relativePermeability(term, perforation->index())/phase->viscosity(term, perforation->index());
        };
        
        injection_perf = std::dynamic_pointer_cast<Injection_Perforation, Perforation>(perforation);
        injection_perf->flow(calculatePeacemanInjection(term,
                                                        main_pressure,
                                                        total_mobility,
                                                        injection_well->injectionPhase(),
                                                        defined_mesh->cell(injection_perf->index()),
                                                        injection_perf->wellIndex(term),
                                                        well->boreholePressure(term),
                                                        well->boreholeDepth()
                                                        ));
    };
    
};

void calculateWellFlow(const int& term, std::shared_ptr<Well>& well){

    double totalFlow = 0;
    
    for(auto perforation = well->begin(); perforation !=well->end(); ++perforation){
        calculatePerforation(term, well, *perforation);
        if((*perforation)->type() == typeid(Production_Perforation).name()){
            
            auto production_perf = std::dynamic_pointer_cast<Production_Perforation, Perforation>(*perforation);
            
            for (auto phase : characterized_phases){
                if(phase->type() != "Gas"){
                    totalFlow += production_perf->flow(phase->index());
                };
            };
            
        }else{
            totalFlow += (*perforation)->totalFlow();
        };
    };

    well->flow(term, totalFlow);
};

#endif /*WELL_FUNCTIONS_H*/
