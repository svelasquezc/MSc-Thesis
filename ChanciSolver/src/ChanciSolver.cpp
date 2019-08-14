#include "NewtonRaphson.h"
#include "Interfluid_Interaction.h"

#include "Producer_Well.h"
#include "Injector_Well.h"


std::vector<std::shared_ptr<Equation_Base>> equations;

std::string timestamp="";
double mytime=0;
double simulationtime = 60;
double timedelta=1;
int Fluid::_count_of_principals=0;
int Fluid::_count_of_fluids=0;
int wells_quantity=0;
int term=0;
int fluids_quantity=0;
int stencil[2] = {-1,1};
int equilibrium_relations_quantity=0;
int cells_number=0;

std::vector<std::shared_ptr<Fluid>> characterized_fluids =
    std::vector<std::shared_ptr<Fluid>>();
std::vector<std::unique_ptr<Equilibrium_Relation>> added_equilibrium_relations =
    std::vector<std::unique_ptr<Equilibrium_Relation>>();

std::vector<std::unique_ptr<Interfluid_Interaction>> added_interfluid_interactions =
    std::vector<std::unique_ptr<Interfluid_Interaction>>();

std::vector<std::shared_ptr<Well>> perforated_wells = std::vector<std::shared_ptr<Well>>();

std::shared_ptr<Mesh> mymesh;
std::shared_ptr<Rock> myrock;

void updateVariables(std::vector<std::shared_ptr<Fluid>>& characterized_fluids, Rock& rock){
    
    for(auto fluid : characterized_fluids){
	fluid->updateProperties(term);
    };

    for(auto equilibrium_relation = added_equilibrium_relations.begin();
	equilibrium_relation!=added_equilibrium_relations.end(); ++equilibrium_relation ){
        
        (equilibrium_relation->get())->updateProperties(term);
            
    };

    rock.updateProperties(term);
};

void calculateInteractions(const int& term, const int& cell_index, Fluid& fluid){

    double capillary_pressure=0;
    
    for(auto interfluid_interaction = added_interfluid_interactions.begin();
	interfluid_interaction!=added_interfluid_interactions.end(); ++interfluid_interaction )
	{
	    auto reference_fluid = interfluid_interaction->get()->referenceFluid();
	    auto wetting_fluid = interfluid_interaction->get()->wettingFluid();
	    auto non_wetting_fluid = interfluid_interaction->get()->nonWettingFluid();
	    
	    if(fluid.index() == reference_fluid->index()){
		
		fluid.relativePermeability(term, cell_index,
					   interfluid_interaction->get()->referenceRelativePermeability(fluid.saturation(term, cell_index)));
		capillary_pressure = interfluid_interaction->get()->capillaryPressure(fluid.saturation(term, cell_index));
	    
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
    

    for(auto interfluid_interaction = added_interfluid_interactions.begin();
	interfluid_interaction!=added_interfluid_interactions.end(); ++interfluid_interaction )
	{
	    auto reference_fluid = interfluid_interaction->get()->referenceFluid();
	    
	    irreducible_saturation = interfluid_interaction->get()->irreducibleSaturation();
	    
	    interpolated_principal_relative_permeability = interfluid_interaction->get()->
		principalRelativePermeability(reference_fluid->saturation(term, cell_index));

	    mobile_saturation = reference_fluid->saturation(term, cell_index) - irreducible_saturation;
	    
	    accumulated_saturation += mobile_saturation;

	    accumulated_principal_relative_permeability += mobile_saturation*interpolated_principal_relative_permeability;

	};
    
    if(accumulated_saturation == 0){
	return 1;
    }else{
	return accumulated_principal_relative_permeability/accumulated_saturation;
    };
    
};

void calculateProperties(const int& term, Cell& cell, Rock& rock){

    const auto cell_index = cell.index();

    double remaining_saturation = 1.0;

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            rock.porosity(term, cell_index, fluid->pressure(term, cell_index));
        }else{
            remaining_saturation = remaining_saturation - fluid->saturation(term, cell_index);

	    calculateInteractions(term, cell_index, *fluid);
        };

    };

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            fluid->saturation(term, cell_index, remaining_saturation);
	    fluid->relativePermeability(term, cell_index, calculateBaker(term, cell_index));
        };

        fluid->volumetricFactor(term, cell_index);
        fluid->viscosity(term, cell_index);

        //Calculate density

        double density_contribution = fluid->standardConditionsDensity();
        
        for(auto equilibrium_relation = added_equilibrium_relations.begin();
            equilibrium_relation!=added_equilibrium_relations.end(); ++equilibrium_relation ){
            
            if((equilibrium_relation->get())->contributorFluid()->index() == fluid->index()){
                
                const auto receiver = (equilibrium_relation->get())->receiverFluid();
                // I suppose this should interpolate the partition coefficient
                // to the contributor fluid pressure
                (equilibrium_relation->get())->partitionCoefficient(term,cell_index);

                density_contribution = density_contribution + 
                    (equilibrium_relation->get())->partitionCoefficient(term,cell_index) *
                    receiver->standardConditionsDensity();
                
            };
        };

        density_contribution = density_contribution/fluid->volumetricFactor(term, cell_index);

        fluid->density(term, cell_index, density_contribution);
        
        fluid->potential(term, cell_index, gravity, cell.depth());
    };        
};

double calculateAccumulation(const int& term, Fluid& fluid, Cell& cell, Rock& rock){

    double past_contribution=0;
    double current_contribution=0;

    const auto cell_index = cell.index();
    
    for(auto equilibrium_relation = added_equilibrium_relations.begin();
        equilibrium_relation!=added_equilibrium_relations.end(); ++equilibrium_relation ){
        
        if((equilibrium_relation->get())->receiverFluid()->index() == fluid.index()){
            
            const auto contributor = (equilibrium_relation->get())->contributorFluid();

            double past_coef = (equilibrium_relation->get())->partitionCoefficient(term-1,cell_index);
            
            past_contribution = past_contribution +
                past_coef * (rock.porosity(term-1,cell_index) * contributor->saturation(term-1,cell_index)
                             / contributor->volumetricFactor(term-1,cell_index));

            double curr_coef = (equilibrium_relation->get())->partitionCoefficient(term,cell_index);
            
            current_contribution = current_contribution +
                curr_coef * (rock.porosity(term,cell_index) * contributor->saturation(term,cell_index)
                             / contributor->volumetricFactor(term,cell_index));
            
        };
    };

    double accumulation = (cell.volume()/timedelta) *
        (((rock.porosity(term,cell_index)*fluid.saturation(term,cell_index)
           /fluid.volumetricFactor(term,cell_index)) + current_contribution)-
         ((rock.porosity(term-1,cell_index)*fluid.saturation(term-1,cell_index)
           /fluid.volumetricFactor(term-1,cell_index)) + past_contribution));
    
    return accumulation;
};

double calculateFlow(const int& term, Fluid& fluid, Mesh& mesh, Cell& cell, Face& face, Rock& rock){
    
    auto harmonicAverage = [](double cell_property, double neighbor_property){
        return 1.0/((1/cell_property) + (1/neighbor_property));
    };

    double flow=0;
    
    int direction = face.orientation();

    const auto neighbor_cell = face.neighbor();
    
    const auto neighbor_index = neighbor_cell->index();
    const auto cell_index = cell.index();

    const auto neighbor_axis = face.neighbor()->numeration3D()[direction];
    const auto axis = cell.numeration3D()[direction];

    double length_delta = mesh.thickness(direction,axis);
    double neighbor_length_delta = mesh.thickness(direction,neighbor_axis);

    double porous_volume = rock.porosity(term, cell_index)*cell.volume();
    double neighbor_porous_volume = rock.porosity(term, neighbor_axis)*neighbor_cell->volume();

    // here is corrected the shape_factor calculation, review the preconceptual schema...
    double shape_factor = face.area()*rock.absolutePermeability(term, cell_index)[direction]/length_delta;
    double neighbor_shape_factor = face.area()*rock.absolutePermeability(term, neighbor_index)[direction]/neighbor_length_delta;

    double face_volumetric_factor = (fluid.volumetricFactor(term, cell_index)*porous_volume +
                                     fluid.volumetricFactor(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_viscosity = (fluid.viscosity(term, cell_index)*porous_volume +
                             fluid.viscosity(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_shape_factor = 2*harmonicAverage(shape_factor, neighbor_shape_factor);

    double face_relative_permeability=1;
    
    //Upwind!!!!
    if(fluid.potential(term, cell_index) >= fluid.potential(term, neighbor_index)){
        face_relative_permeability = fluid.relativePermeability(term, cell_index);
    }else{
        face_relative_permeability = fluid.relativePermeability(term, neighbor_index);
    };

    double face_transmissivity = face_shape_factor * face_relative_permeability /
        (face_volumetric_factor * face_viscosity);

    flow = flow + face_transmissivity*(fluid.potential(term, cell_index) - fluid.potential(term, neighbor_index));

    for(auto equilibrium_relation = added_equilibrium_relations.begin();
        equilibrium_relation!=added_equilibrium_relations.end(); ++equilibrium_relation ){
        
        if((equilibrium_relation->get())->receiverFluid()->index() == fluid.index()){

            const auto contributor = (equilibrium_relation->get())->contributorFluid();
            
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
                ((equilibrium_relation->get())->partitionCoefficient(term, cell_index)*porous_volume +
                 (equilibrium_relation->get())->partitionCoefficient(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);
            
            double face_transmissivity = face_shape_factor * face_relative_permeability /
                (face_volumetric_factor * face_viscosity);

            flow = flow + face_partition_coefficient*face_transmissivity*
                (contributor->potential(term, cell_index) - contributor->potential(term, neighbor_index));
        };
    };

    return flow;
};

using BlackOilNewton = NewtonRaphson<decltype(calculateProperties), decltype(calculateFlow), decltype(calculateAccumulation)>;

BlackOilNewton my_newton(calculateProperties,calculateFlow,calculateAccumulation);


//Change Event Name
void FluidPressureVaries(std::string& _timestamp){
    if(_timestamp == "stop"){
        updateVariables(characterized_fluids, *myrock);
    }
};

void timePasses(std::string& _timestamp, int& _term, double& _mytime, double& _timedelta, double& _simulationtime){
    if(_timestamp == "continue" && _mytime<=_simulationtime){
        std::cout << "Simulating interval [" << _mytime << " - " << _mytime + _timedelta << "]" << std::endl;
        
        _mytime    +=_timedelta;
        _timestamp = "stop";
        ++_term;      
    }
};

//Rock
void launchTriggers(){
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    mymesh->appear(timestamp,stencil);
    
};

void launchGeomodeler(){
    int option;
    int _dimension;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Define Mesh";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        mymesh = std::make_shared<Mesh>(Mesh());
        mymesh->defineMesh();
        cells_number = mymesh->getCellTotal();
        break;
    default:
        break;
    }
}

void launchPetrophysicalEngineer(){
    int option;    
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Rock";
    std::cout << "2. Add Interfluid Interaction";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        myrock = std::make_shared<Rock>(Rock());
        myrock->characterize(cells_number);
        break;
    case 2:
        if(fluids_quantity >= 2){
            added_interfluid_interactions.push_back(std::make_unique<Interfluid_Interaction>());
            (--added_interfluid_interactions.end())->get()->add(fluids_quantity,characterized_fluids);
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
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Fluid" << std::endl;
    std::cout << "2. Add Equilibrium Relation";
    
    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        characterized_fluid = std::make_shared<Fluid>(Fluid());
        characterized_fluids.push_back(characterized_fluid);
        (--(characterized_fluids.end()))->get()->characterize(cells_number);
        equations.push_back(characterized_fluid);
        ++fluids_quantity;
        break;
    case 2:
        if(fluids_quantity >= 2){
            added_equilibrium_relations.push_back(std::make_unique<Equilibrium_Relation>());
            (--(added_equilibrium_relations.end()))->get()->add(fluids_quantity,characterized_fluids);
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

void launchReservoirEngineer(){
    int option;
    std::string type;
    
    std::shared_ptr<Well> well;
    
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Perforate Well" << std::endl;

        Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        if(fluids_quantity >= 1){
            Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
            if(type == "Producer"){
                well = std::make_shared<Producer_Well>();
            }else{
                well = std::make_shared<Injector_Well>();
            }

            
        }else{
            std::cout << "It is not possible to perforate wells with no fluid characterized."
                      << std::endl;
        };
        break;
    default:
        break;
    }
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
            launchReservoirEngineer();
            break;
        case -1:
            run=true;
            break;
        }
        launchTriggers();
    };
}

int main(){
    launchMenu();
    timestamp="continue";
    
    while(mytime<simulationtime){
        launchTriggers();
    };
    
    return 0;
};
