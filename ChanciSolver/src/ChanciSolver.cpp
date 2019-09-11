#include <ctime>

#include "Global.h"

#include "NewtonRaphson.h"
#include "Equilibrium_Relation.h"
#include "Interfluid_Interaction.h"


int Fluid::_count_of_principals=0;
int Fluid::_count_of_fluids=0;

std::vector<std::shared_ptr<Equation_Base>> equations =
    std::vector<std::shared_ptr<Equation_Base>>();

std::vector<std::shared_ptr<Fluid>> characterized_fluids =
    std::vector<std::shared_ptr<Fluid>>();

std::vector<std::unique_ptr<Equilibrium_Relation>> added_equilibrium_relations =
    std::vector<std::unique_ptr<Equilibrium_Relation>>();

std::vector<std::unique_ptr<Interfluid_Interaction>> added_interfluid_interactions =
    std::vector<std::unique_ptr<Interfluid_Interaction>>();

std::vector<std::shared_ptr<Well>> perforated_wells =
    std::vector<std::shared_ptr<Well>>();

std::unique_ptr<Mesh> mymesh;
std::unique_ptr<Rock> myrock;

void updateVariables(const int& term, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, Rock& rock){
    
    for(auto fluid : characterized_fluids){
	fluid->updateProperties(term);
    };

    for(auto& equilibrium_relation : added_equilibrium_relations ){
        equilibrium_relation->updateProperties(term);            
    };

    for(auto well : perforated_wells){
        well->updateProperties(term);
    };

    rock.updateProperties(term);
};

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

void calculateInteractions(const int& term, const int& cell_index, Fluid& fluid){

    double capillary_pressure=0;
    
    for(auto& interfluid_interaction : added_interfluid_interactions )
        {
            auto reference_fluid = interfluid_interaction->referenceFluid();
            auto wetting_fluid = interfluid_interaction->wettingFluid();
            auto non_wetting_fluid = interfluid_interaction->nonWettingFluid();
	    
            if(fluid.index() == reference_fluid->index()){
		
                fluid.relativePermeability(term, cell_index,
                                           interfluid_interaction->referenceRelativePermeability(fluid.saturation(term, cell_index)));
                capillary_pressure = interfluid_interaction->capillaryPressure(fluid.saturation(term, cell_index));
	    
                if(fluid.index() == wetting_fluid->index()){
                    fluid.pressure(term, cell_index,
                                   non_wetting_fluid->pressure(term, cell_index) - capillary_pressure);
                }else{
                    fluid.pressure(term, cell_index,
                                   wetting_fluid->pressure(term, cell_index) + capillary_pressure);
                };
            };
        };
};

double calculateBaker(const int& term, const int& cell_index){

    double accumulated_saturation = 0;
    double accumulated_principal_relative_permeability=0;
    double irreducible_saturation = 0;
    double mobile_saturation=0;
    double interpolated_principal_relative_permeability=0;
    double maximum_principal_relative_permeability=0;
    

    for(auto& interfluid_interaction : added_interfluid_interactions )
        {
            auto reference_fluid = interfluid_interaction->referenceFluid();
	    
            irreducible_saturation = interfluid_interaction->irreducibleSaturation();
	    
            interpolated_principal_relative_permeability = interfluid_interaction->
                principalRelativePermeability(reference_fluid->saturation(term, cell_index));

            mobile_saturation = reference_fluid->saturation(term, cell_index) - irreducible_saturation;
	    
            accumulated_saturation += mobile_saturation;

            accumulated_principal_relative_permeability += mobile_saturation*interpolated_principal_relative_permeability;

            maximum_principal_relative_permeability = std::max(maximum_principal_relative_permeability, interfluid_interaction->maximumPrincipalRelativePermeability());
            
        };
    
    if(accumulated_saturation == 0){
        return maximum_principal_relative_permeability;
    }else{
        return accumulated_principal_relative_permeability/accumulated_saturation;
    };
    
};

void calculateProperties(const int& term, const std::shared_ptr<Cell>& cell, Rock& rock){

    const auto cell_index = cell->index();

    double remaining_saturation = 1.0;

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            rock.porosity(term, cell_index, fluid->pressure(term, cell_index));
            //rock.porosity(term, cell_index, fluid->pressure(term, cell_index));
        }else{
            remaining_saturation = remaining_saturation - fluid->saturation(term, cell_index);

            calculateInteractions(term, cell_index, *fluid);
        };

    };

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            fluid->saturation(term, cell_index, remaining_saturation);
            if(Global::fluids_quantity>=2){
                fluid->relativePermeability(term, cell_index, calculateBaker(term, cell_index));
            }else{
                fluid->relativePermeability(term, cell_index, 1.0);
            };
        };

        fluid->volumetricFactor(term, cell_index, fluid->pressure(term, cell_index));
        fluid->viscosity(term, cell_index, fluid->pressure(term, cell_index));
        //Calculate density

        double density_contribution = fluid->standardConditionsDensity();
        
        for(auto& equilibrium_relation : added_equilibrium_relations ){
            
            if(equilibrium_relation->contributorFluid()->index() == fluid->index()){
                
                const auto receiver = equilibrium_relation->receiverFluid();
                // I suppose this should interpolate the partition coefficient
                // to the contributor fluid pressure
                equilibrium_relation->partitionCoefficient(term,cell_index);

                density_contribution = density_contribution + 
                    equilibrium_relation->partitionCoefficient(term,cell_index) *
                    receiver->standardConditionsDensity();
                
            };
        };

        density_contribution = density_contribution/fluid->volumetricFactor(term, cell_index);

        fluid->density(term, cell_index, density_contribution);
        
        fluid->potential(term, cell_index, gravity, cell->depth());
    };        
};

double calculateAccumulation(const int& term, Fluid& fluid, const std::shared_ptr<Cell>& cell, Rock& rock){

    double past_contribution=0;
    double current_contribution=0;

    //std::string N="N";
    //std::string K="K";

    const int cell_index = cell->index();
    
    for(auto& equilibrium_relation : added_equilibrium_relations ){
        
        if(equilibrium_relation->receiverFluid()->index() == fluid.index()){
            
            const auto contributor = equilibrium_relation->contributorFluid();

            double past_coef = equilibrium_relation->partitionCoefficient(term-1,cell_index);
            
            past_contribution = past_contribution +
                past_coef * (rock.porosity(term-1,cell_index) * contributor->saturation(term-1,cell_index)
                             / contributor->volumetricFactor(term-1,cell_index));

            double curr_coef = equilibrium_relation->partitionCoefficient(term,cell_index);
            
            current_contribution = current_contribution +
                curr_coef * (rock.porosity(term,cell_index) * contributor->saturation(term,cell_index)
                             / contributor->volumetricFactor(term,cell_index));
            
        };
    };

    double accumulation = (cell->volume()/Global::timedelta) *
        (((rock.porosity(term,cell_index)*fluid.saturation(term,cell_index)
           /fluid.volumetricFactor(term,cell_index)) + current_contribution)-
         ((rock.porosity(term-1,cell_index)*fluid.saturation(term-1,cell_index)
           /fluid.volumetricFactor(term-1,cell_index)) + past_contribution));
    
    return accumulation;
};

double calculateFlow(const int& term, Fluid& fluid, const Mesh& mesh, const std::shared_ptr<Cell>& cell, const std::shared_ptr<Face>& face, Rock& rock){
    
    auto harmonicAverage = [](double cell_property, double neighbor_property){
        return 1.0/((1.0/cell_property) + (1.0/neighbor_property));
    };

    double flow=0;
    
    int direction = face->orientation();

    auto neighbor_cell = face->neighbor().lock();
    
    const auto neighbor_index = neighbor_cell->index();
    const auto cell_index = cell->index();

    const auto neighbor_axis = neighbor_cell->numeration3D()[direction];
    const auto axis = cell->numeration3D()[direction];

    double length_delta = mesh.thickness(direction,axis);
    double neighbor_length_delta = mesh.thickness(direction,neighbor_axis);

    double porous_volume = rock.porosity(term, cell_index)*cell->volume();
    double neighbor_porous_volume = rock.porosity(term, neighbor_axis)*neighbor_cell->volume();

    // here is corrected the shape_factor calculation, review the preconceptual schema...
    double shape_factor = face->area()*rock.absolutePermeability(term, cell_index)[direction]/length_delta;
    double neighbor_shape_factor = face->area()*rock.absolutePermeability(term, neighbor_index)[direction]/neighbor_length_delta;

    double face_volumetric_factor = (fluid.volumetricFactor(term, cell_index)*porous_volume +
                                     fluid.volumetricFactor(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_viscosity = (fluid.viscosity(term, cell_index)*porous_volume +
                             fluid.viscosity(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_shape_factor = 2.0*harmonicAverage(shape_factor, neighbor_shape_factor);

    double face_relative_permeability=1;
    
    //Upwind!!!!
    if(fluid.potential(term, cell_index) >= fluid.potential(term, neighbor_index)){
        face_relative_permeability = fluid.relativePermeability(term, cell_index);
    }else{
        face_relative_permeability = fluid.relativePermeability(term, neighbor_index);
    };

    double face_transmissivity = face_shape_factor * face_relative_permeability /
        (face_volumetric_factor * face_viscosity);

    flow = flow + face_transmissivity*(fluid.potential(term, neighbor_index) - fluid.potential(term, cell_index));

    for(auto& equilibrium_relation : added_equilibrium_relations ){
        
        if(equilibrium_relation->receiverFluid()->index() == fluid.index()){

            const auto contributor = equilibrium_relation->contributorFluid();
            
            double face_volumetric_factor = (contributor->volumetricFactor(term, cell_index)*porous_volume +
                                             contributor->volumetricFactor(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);

            double face_viscosity = (contributor->viscosity(term, cell_index)*porous_volume +
                                     contributor->viscosity(term, neighbor_index)*neighbor_porous_volume) /
                (porous_volume + neighbor_porous_volume);

            //Upwind!!!!
            if(contributor->potential(term, cell_index) >= contributor->potential(term, neighbor_index)){
                face_relative_permeability = contributor->relativePermeability(term, cell_index);
            }else{
                face_relative_permeability = contributor->relativePermeability(term, neighbor_index);
            };

            double face_partition_coefficient =
                (equilibrium_relation->partitionCoefficient(term, cell_index)*porous_volume +
                 equilibrium_relation->partitionCoefficient(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);
            
            double face_transmissivity = face_shape_factor * face_relative_permeability /
                (face_volumetric_factor * face_viscosity);

            flow = flow + face_partition_coefficient*face_transmissivity*
                (contributor->potential(term, neighbor_index) - contributor->potential(term, cell_index));
        };
    };

    neighbor_cell.reset();
    
    return flow;
};

using BlackOilNewton = NewtonRaphson<decltype(calculateProperties), decltype(calculateFlow), decltype(calculateAccumulation),decltype(calculatePerforation),decltype(calculateWellFlow), decltype(estimateWellPressure)>;

std::unique_ptr<BlackOilNewton> my_newton;// = BlackOilNewton(0,0,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

//Change Event Name
void FluidPressureVaries(std::string& timestamp){
    if(timestamp == "stop"){

        updateVariables(1,characterized_fluids, *myrock);

        my_newton->iterate(Global::term, *mymesh, perforated_wells, equations, *myrock);

        updateVariables(0,characterized_fluids, *myrock);
        
        timestamp = "continue";

        double next_time = Global::mytime + Global::timedelta;

        for(auto well : perforated_wells){
            if(next_time >= well->operativeCondition()->nextChange()){
                ++Global::changing_wells;
                timestamp == "change";
                well->operativeStatus(2); //Pending Change
            };
        };
    }
};

void timePasses(std::string& timestamp, int& term, double& mytime, double& timedelta, double& simulationtime){
    if(timestamp == "continue" && mytime<=simulationtime){
        int max_number_of_well_non_zeros=0;
        int well_equations=0;
        int total_equations;
        int max_non_zeros;
        if(mytime == 0){
            for(auto cell = mymesh->begin(); cell != mymesh->end(); ++cell){
                calculateProperties(0, *cell, *myrock);
            };

            for(auto well : perforated_wells){
                
                if(well->operativeStatus() == 1){ //Just Changed

                    if(well->operativeCondition()->type() == "FLOW"){
                        estimateWellPressure(0, well);
                        max_number_of_well_non_zeros += well->numberOfPerforates()*2+1;
                        ++well_equations;
                    }else{
                        calculateWellFlow(0, well);
                    };
                    
                    well->operativeStatus(0); // Stable
                };
            };

            total_equations = Global::fluids_quantity*Global::cells_number + well_equations;

            max_non_zeros = Global::fluids_quantity*Global::fluids_quantity*Global::cells_number + max_number_of_well_non_zeros;
            
            my_newton = std::make_unique<BlackOilNewton>(total_equations,max_non_zeros,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);
            
            ++term;
            
        };

        if(Global::changing_wells > 0){
            
            max_number_of_well_non_zeros=0;
            well_equations=0;

            for(auto well : perforated_wells){
                
                if(well->operativeCondition()->type() == "FLOW"){
                    max_number_of_well_non_zeros += well->numberOfPerforates()*2+1;
                    ++well_equations;
                };                    

            };
            
            total_equations = Global::fluids_quantity*Global::cells_number + well_equations;
            max_non_zeros = Global::fluids_quantity*Global::fluids_quantity*Global::cells_number + max_number_of_well_non_zeros;
            my_newton.reset();
            my_newton = std::make_unique<BlackOilNewton>(total_equations,max_non_zeros,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

            Global::changing_wells = 0;
        };
        
        std::cout << "Simulating interval [" << mytime << " - " << mytime + timedelta << "]" << std::endl;

        for (auto cell = mymesh->begin(); cell!=mymesh->end(); ++cell){
            std::cout << characterized_fluids[0]->pressure(0, (*cell)->index())<< " ";
        };
        std::cout << std::endl;
        
        mytime    +=timedelta;
        timestamp = "stop";
    };
};

//Rock

void launchGeomodeler(){
    int option;
    int _dimension;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Define Mesh";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        mymesh = std::make_unique<Mesh>();
        mymesh->define();
        Global::cells_number = mymesh->getCellTotal();
        break;
    default:
        break;
    }
}

void launchPetrophysicalEngineer(){
    int option;
    std::unique_ptr<Interfluid_Interaction> added_interfluid_interaction;
    
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Rock" << std::endl;;
    std::cout << "2. Add Interfluid Interaction";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        myrock = std::make_unique<Rock>();
        myrock->characterize(Global::cells_number);
        break;
    case 2:
        if(Global::fluids_quantity >= 2){
            
            added_interfluid_interaction = std::make_unique<Interfluid_Interaction>();
            added_interfluid_interaction->add(characterized_fluids);
            
            added_interfluid_interactions.push_back(std::move(added_interfluid_interaction));
            
        }else{
            std::cout << "It is not possible to add an Interfluid interaction with only one fluid characterized."
                      << std::endl;
        }
        break;
    default:
        
        break;
    }
};

void launchFluidsEngineer(){
    int option;
    int _dimension;
    std::shared_ptr<Fluid> characterized_fluid;
    std::unique_ptr<Equilibrium_Relation> added_equilibrium_relation;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Fluid" << std::endl;
    std::cout << "2. Add Equilibrium Relation";
    
    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        characterized_fluid = std::make_shared<Fluid>(Global::fluids_quantity);
        characterized_fluid->characterize(Global::cells_number);
        characterized_fluids.push_back(characterized_fluid);
        equations.push_back(characterized_fluid);
        ++Global::fluids_quantity;
        break;
    case 2:
        if(Global::fluids_quantity >= 2){
            added_equilibrium_relation = std::make_unique<Equilibrium_Relation>(Global::equilibrium_relations_quantity);
            added_equilibrium_relation->add(Global::cells_number,characterized_fluids);
            added_equilibrium_relations.push_back(std::move(added_equilibrium_relation));
            ++Global::equilibrium_relations_quantity;
            break;
        }else{
            std::cout << "It is not possible to add an Equilibrium relation with only one fluid characterized."
                      << std::endl;
        }
        break;
    default:
        break;
    }

};

void launchReservoirEngineer(std::string& timestamp){
    int option;
    std::string type;
    
    std::shared_ptr<Well> well;

    int index=0;
    
    std::cout << "Select your action" << std::endl;
    if(timestamp =="")std::cout << "1. Perforate Well" << std::endl;
    std::cout << "2. Establish Operative Condition ";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    if(timestamp =="change"){
        while(option != 2){
            Value_Reader::myRead(std::string("Please insert a valid option"), option, std::string("Please insert a valid option"));
        };
    };
    
    switch(option){
    case 1:
        if(timestamp ==""){
            if(Global::fluids_quantity >= 1){
                Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                while(type != "Producer" && type != "Injector"){
                    std::cout << "Only Injector or Producer wells";
                    Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                };
                ++Global::wells_quantity;
                if(type == "Producer"){
                    well = std::make_shared<Producer_Well>(Producer_Well(Global::wells_quantity));
                }else if(type == "Injector"){
                    well = std::make_shared<Injector_Well>(Injector_Well(Global::wells_quantity));
                }else{
                    std::cout << "Only Injector or Producer wells";
                };
            
                well->perforate(*mymesh, characterized_fluids,type);
                perforated_wells.push_back(well);
                equations.push_back(well);
            
            }else{
                std::cout << "It is not possible to perforate wells with no fluid characterized."
                          << std::endl;
            };
        };
        break;
    case 2:
        std::cout << "Please select the index of the well:"
                  << std::endl;
        for (auto well : perforated_wells){
            std::cout << well->index() <<". " <<well->type()<<std::endl;
        };
        while(index<1 && index>=Global::wells_quantity){
            Value_Reader::myRead(std::string("Please select an index between the range "), index, std::string("Please insert a valid input"));
        };

        perforated_wells[index-1]->establish(Global::term, Global::timestamp);
        
    default:
        break;
    }
};

void reEstablishOperativeConditions(std::ifstream& file_reader, const int& term, std::string& timestamp){
    
    std::string object;
    
    if(timestamp == "change"){
        for(auto well : perforated_wells){
            if(well->operativeStatus() == 2){ //Pending Change
                file_reader >> object;
                std::transform(object.begin(), object.end(),object.begin(), ::toupper);
                if(object == "WELL"){
                    int index;
                    file_reader >> index;
                    if(well->index() == index){
                        well->establishFromFile(file_reader, term, timestamp);
                    };
                };
            };
        };
        timestamp == "continue";
    };
};

void launchTriggers(){
    using namespace Global;
    mymesh->appear(timestamp,stencil);
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    FluidPressureVaries(timestamp);
    if(timestamp == "change"){
        launchReservoirEngineer(timestamp);
    }
};

void launchTriggers(std::ifstream& file_reader){
    using namespace Global;
    mymesh->appear(timestamp,stencil);
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    FluidPressureVaries(timestamp);
    reEstablishOperativeConditions(file_reader, term, timestamp);
};

void launchMenu(){
    int option;
    bool run=false;
    while(!run){
        std::cout << "Select your role or -1 for running simulation" << std::endl;
        std::cout << "1. Geomodeler"<< std::endl << "2. Petrophysical Engineer" << std::endl;
        std::cout << "3. Fluids Engineer" << std::endl << "4. Reservoir Engineer";
        
        Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid role"));
        
        switch(option){
        case 1:
            launchGeomodeler();
            break;
        case 2:
            launchPetrophysicalEngineer();
            break;
        case 3:
            launchFluidsEngineer();
            break;
        case 4:
            launchReservoirEngineer(Global::timestamp);
            break;
        case -1:
            run=true;
            break;
        }
        launchTriggers();
    };
};


void launchFromFile(std::ifstream& file_reader){
    
    std::string object;
    std::string type;
    
    std::shared_ptr<Fluid> characterized_fluid;
    std::unique_ptr<Equilibrium_Relation> added_equilibrium_relation;
    std::unique_ptr<Interfluid_Interaction> added_interfluid_interaction;
    std::shared_ptr<Well> well;
    
    if(file_reader.is_open()){
        while(file_reader>>object){
            
            std::transform(object.begin(), object.end(),object.begin(), ::toupper);
            
            if(object == "MESH"){
                
                mymesh = std::make_unique<Mesh>();
                mymesh->defineFromFile(file_reader);
                Global::cells_number = mymesh->getCellTotal();
                
            }else if(object == "ROCK"){
                
                myrock = std::make_unique<Rock>();
                myrock->characterizeFromFile(file_reader, Global::cells_number);
            
            }else if(object == "FLUIDS"){
                
                file_reader>>Global::fluids_quantity;
                
                for(int fluid=0; fluid<Global::fluids_quantity;++fluid){
                    
                    characterized_fluid = std::make_shared<Fluid>(fluid);
                    characterized_fluid->characterizeFromFile(file_reader, Global::cells_number);
                    characterized_fluids.push_back(characterized_fluid);
                    equations.push_back(characterized_fluid);
                    
                };
                
            }else if(object == "EQUILIBRIUM_RELATIONS"){
                
                if(Global::fluids_quantity >= 2){
                    file_reader>>Global::equilibrium_relations_quantity;
                    for(int equilibrium_relation=0; equilibrium_relation<Global::equilibrium_relations_quantity; ++equilibrium_relation){
                        added_equilibrium_relation = std::make_unique<Equilibrium_Relation>(equilibrium_relation);
                        added_equilibrium_relation->addFromFile(file_reader, Global::cells_number, characterized_fluids);
                        added_equilibrium_relations.push_back(std::move(added_equilibrium_relation));
                    };
                }else{
                    std::cout << "There is only one fluid, EQUILIBRIUM_RELATIONS will be omitted\n";
                };
                
            }else if(object == "INTERFLUID_INTERACTIONS"){
                
                if(Global::fluids_quantity >= 2){
                    file_reader>>Global::interfluid_interactions_quantity;
                    for(int interfluid_interaction=0; interfluid_interaction<Global::interfluid_interactions_quantity; ++interfluid_interaction){
                        added_interfluid_interaction = std::make_unique<Interfluid_Interaction>(interfluid_interaction);
                        added_interfluid_interaction->addFromFile(file_reader, characterized_fluids);
                        added_interfluid_interactions.push_back(std::move(added_interfluid_interaction));
                    };
                }else{
                    std::cout << "There is only one fluid, INTERFLUID_INTERACTIONS will be omitted\n";
                };
                
            }else if(object == "TIME_DELTA"){
                
                file_reader>>Global::timedelta;
                
            }else if(object == "SIMULATION_TIME"){
                
                file_reader>>Global::simulationtime;
                
            }else if(object == "WELLS"){
                if(Global::fluids_quantity >= 1){
                    file_reader>>Global::wells_quantity;
                    for(int well_index=0; well_index<Global::wells_quantity;++well_index){
                        file_reader>>type;
                        std::transform(type.begin(), type.end(), type.begin(), ::toupper);
                        if(type == "PRODUCER"){
                            well = std::make_shared<Producer_Well>(Producer_Well(well_index));
                        }else if(type == "INJECTOR"){
                            well = std::make_shared<Injector_Well>(Injector_Well(well_index));
                        }else{
                            
                            std::ostringstream well_err = std::ostringstream();
                            
                            well_err<<"The type selected for well "<< well_index+1 <<": ("<< type <<") is invalid, only Producer or Injector Wells"<<std::endl;
                            
                            throw std::invalid_argument(well_err.str());
                        };
                        well->perforateFromFile(file_reader,*mymesh, characterized_fluids,type);
                        perforated_wells.push_back(well);
                        equations.push_back(well);

                    };
                }else{
                    throw std::domain_error("There are no fluids characterized or WELLS Keyword appears before FLUIDS");
                };
                
            }else if(object == "OPERATIVE_CONDITIONS"){
                if(Global::wells_quantity >= 1){
                    for(auto well : perforated_wells){
                        file_reader >> object;
                        std::transform(object.begin(), object.end(),object.begin(), ::toupper);
                        if(object == "WELL"){
                            int index;
                            file_reader >> index;
                            if(well->index() == index - 1){
                                well->establishFromFile(file_reader, Global::term, Global::timestamp);
                            };
                        };
                    };
                    break;
                }else{
                    throw std::domain_error("There are no wells perforated or OPERATIVE_CONDITIONS appears before WELLS");
                };
            };
            launchTriggers(file_reader);
        };
        file_reader.close();
    }else{
        std::cout << "No success opening the file"<<std::endl;
    };
    
};

int main(int argc, char *argv[]){
    std::time_t tstart;
    std::time_t tend;
    if(argc <= 1){
        launchMenu();
        
        Global::timestamp="continue";
        tstart = std::time(0);
        while(Global::mytime<Global::simulationtime){
            launchTriggers();
        };
        tend = std::time(0);
        
        std::cout << "Elapsed time: "<< std::difftime(tend, tstart) <<" seconds(s)."<<std::endl;
        
    }else{
        std::ifstream file_reader(argv[1], std::ifstream::in);
        launchFromFile(file_reader);
        Global::timestamp="continue";
        tstart = std::time(0);
        while(Global::mytime<Global::simulationtime){
            launchTriggers(file_reader);
        };
        tend = std::time(0);

        std::cout << "Elapsed time: "<< std::difftime(tend, tstart) <<" seconds(s)."<<std::endl;
    };
    
    return 0;
};
