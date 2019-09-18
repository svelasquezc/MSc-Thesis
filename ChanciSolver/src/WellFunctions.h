#ifndef WELL_FUNCTIONS_H
#define WELL_FUNCTIONS_H

#include "Producer_Well.h"
#include "Injector_Well.h"

void calculateGeometry(const int& term, std::shared_ptr<Well>& well, std::shared_ptr<Perforate>& perforation){
    
    auto cell = mymesh->cell(perforation->index());

    auto numeration = cell->numeration3D();

    auto absolute_permeability= myrock->absolutePermeability(term,cell->index());
    
    int axis_x = numeration[0];
    int axis_y = numeration[1];
    int axis_z = numeration[2];

    perforation->calculateEquivalentRadius(term, absolute_permeability[1],
                                           absolute_permeability[0],
                                           mymesh->thickness(1,axis_y),
                                           mymesh->thickness(0,axis_x)
                                           );
    
    perforation->calculateWellIndex       (term,
                                           mymesh->thickness(2,axis_z),
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
    
    
    std::shared_ptr<Injector_Perforate> injector_perf;
    std::shared_ptr<Producer_Perforate> producer_perf;

    for(auto perforation = well->begin(); perforation !=well->end(); ++perforation){    
        calculateGeometry(term, well, *perforation);
    };

    std::shared_ptr<Injector_Well> injector_well;
    std::shared_ptr<Producer_Well> producer_well;

    std::shared_ptr<Fluid> main_fluid;

    for(auto fluid : characterized_fluids){
        if(fluid->principal()){
            main_fluid = fluid;
        };
    }
    
    if(well->type() == typeid(Producer_Well).name()){
        producer_well = std::dynamic_pointer_cast<Producer_Well, Well>(well);
        for(auto perforation = producer_well->begin(); perforation !=producer_well->end(); ++perforation){
            producer_perf = std::dynamic_pointer_cast<Producer_Perforate, Perforate>(*perforation);

            auto cell = mymesh->cell(producer_perf->index());
            
            for(auto fluid : characterized_fluids){
                
                if(fluid->type() != "Gas"){
                    perforation_mobility = producer_perf->wellIndex(term)*
                        fluid->relativePermeability(term, cell->index())/(fluid->volumetricFactor(term, cell->index())*fluid->viscosity(term, cell->index()));

                    accumulated_mobility += perforation_mobility;

                    hydrostatic_head = fluid->density(term, cell->index())*gravity*(producer_well->boreholeDepth()-cell->depth());

                    estimated_flow = perforation_mobility*(main_fluid->pressure(term, cell->index()) + hydrostatic_head);
                    accumulated_flow += estimated_flow;
                    
                };
                
            };
            
        };
        
    }else{
        
        injector_well = std::dynamic_pointer_cast<Injector_Well, Well>(well);
        auto injection_fluid = injector_well->injectionFluid();
        
        for(auto perforation = injector_well->begin(); perforation !=injector_well->end(); ++perforation){
            injector_perf = std::dynamic_pointer_cast<Injector_Perforate, Perforate>(*perforation);

            auto cell = mymesh->cell(injector_perf->index());

            double total_mobility = 0;

            std::shared_ptr<Fluid> main_fluid;
        
            for(auto fluid : characterized_fluids){
                if(fluid->principal()){
                    main_fluid = fluid;
                };
                total_mobility += fluid->relativePermeability(term, cell->index())/fluid->viscosity(term, cell->index());
            };
            
            perforation_mobility = injector_perf->wellIndex(term)*total_mobility
                /injection_fluid->volumetricFactor(term, cell->index());

            accumulated_mobility += perforation_mobility;

            hydrostatic_head = injection_fluid->density(term, cell->index())*gravity*(injector_well->boreholeDepth()-cell->depth());

            estimated_flow = perforation_mobility*(main_fluid->pressure(term, cell->index()) + hydrostatic_head);
            accumulated_flow += estimated_flow;
                    
        };
        
    };

    well->boreholePressure(term, (accumulated_flow + well->flow(term))/accumulated_mobility);
};

double calculatePeacemanProducer(const int& term, const double main_pressure, const std::shared_ptr<Fluid>& fluid, const std::shared_ptr<Cell>& cell, const double well_index, const double borehole_pressure, const double borehole_depth){
    
    auto cell_index = cell->index();
    double peaceman_flow = 0.0;
    
    peaceman_flow = (well_index * fluid->relativePermeability(term, cell_index) / ((fluid->volumetricFactor(term, cell_index)*fluid->viscosity(term, cell_index)))) *
        (borehole_pressure - main_pressure - (fluid->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));
    
    //std::cout << "\n mira::wi  " <<    well_index ;
    //std::cout << "\n mira::kr  " << fluid->relativePermeability(term, cell_index);
    //std::cout << "\n mira::bol  " << fluid->volumetricFactor(term, cell_index);
    //std::cout << "\n mira::vis  " << fluid->viscosity(term, cell_index);
    //std::cout << "\n mira::bhp  " << borehole_pressure;
    //std::cout << "\n mira::bhp  " <<    main_pressure;
    //std::cout << "\n mira::Pres  " << (borehole_pressure - main_pressure - (fluid->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));

    
    //std::cout <<"\n mira:" << borehole_depth - cell->depth() << "\n";
    
};

    double calculatePeacemanInjector(const int& term, const double main_pressure, const double total_mobility, const std::shared_ptr<Fluid>& fluid, const std::shared_ptr<Cell>& cell, const double well_index, const double borehole_pressure, const double borehole_depth){
    
    auto cell_index = cell->index();
    double peaceman_flow = 0;

    peaceman_flow = (well_index * total_mobility / fluid->volumetricFactor(term, cell_index)) *
        (borehole_pressure - main_pressure - (fluid->density(term, cell_index)*gravity*(borehole_depth - cell->depth())));
    
};

void calculatePerforation(const int& term, std::shared_ptr<Well>& well, std::shared_ptr<Perforate>& perforation){

    std::shared_ptr<Injector_Perforate> injector_perf;
    std::shared_ptr<Producer_Perforate> producer_perf;

    
    std::shared_ptr<Injector_Well> injector_well;
    std::shared_ptr<Producer_Well> producer_well;

    if(well->type() == typeid(Producer_Well).name()){
        producer_well = std::dynamic_pointer_cast<Producer_Well, Well>(well);
    }else{
        injector_well = std::dynamic_pointer_cast<Injector_Well, Well>(well);
    };

    double main_pressure;

    for(auto fluid : characterized_fluids){
        if(fluid->principal()){
            main_pressure = fluid->pressure(term, perforation->index());
        }
    }

    calculateGeometry(term, well, perforation);
    
    if(perforation->type() == typeid(Producer_Perforate).name()){
        
        producer_perf = std::dynamic_pointer_cast<Producer_Perforate, Perforate>(perforation);
        
        for(auto fluid : characterized_fluids){
            producer_perf->flow(fluid->index(),
                                calculatePeacemanProducer(term,
                                                          main_pressure,
                                                          fluid,
                                                          mymesh->cell(producer_perf->index()),
                                                          producer_perf->wellIndex(term),
                                                          well->boreholePressure(term),
                                                          well->boreholeDepth()
                                                          ));
        };
        
    }else{
        
        double total_mobility = 0;
        
        for(auto fluid : characterized_fluids){
            total_mobility += fluid->relativePermeability(term, perforation->index())/fluid->viscosity(term, perforation->index());
        };
        
        injector_perf = std::dynamic_pointer_cast<Injector_Perforate, Perforate>(perforation);
        injector_perf->flow(calculatePeacemanInjector(term,
                                                      main_pressure,
                                                      total_mobility,
                                                      injector_well->injectionFluid(),
                                                      mymesh->cell(injector_perf->index()),
                                                      injector_perf->wellIndex(term),
                                                      well->boreholePressure(term),
                                                      well->boreholeDepth()
                                                      ));
    };
    
};

void calculateWellFlow(const int& term, std::shared_ptr<Well>& well){

    double totalFlow = 0;
    
    for(auto perforation = well->begin(); perforation !=well->end(); ++perforation){
        calculatePerforation(term, well, *perforation);
        if((*perforation)->type() == typeid(Producer_Perforate).name()){
            
            auto producer_perf = std::dynamic_pointer_cast<Producer_Perforate, Perforate>(*perforation);
            
            for (auto fluid : characterized_fluids){
                if(fluid->type() != "Gas"){
                    totalFlow += producer_perf->flow(fluid->index());
                };
            };
            
        }else{
            totalFlow += (*perforation)->totalFlow();
        };
    };

    well->flow(term, totalFlow);
};

#endif /*WELL_FUNCTIONS_H*/
