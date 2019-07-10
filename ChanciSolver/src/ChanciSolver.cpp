
#include <cmath>

#include "NewtonRaphson.h"

std::string timestamp="";
double mytime=0;
double simulationtime = 60;
double timedelta=1;
int term=0;
int fluids_quantity=0;
int stencil[2] = {-1,1};
int equilibrium_relations_quantity=0;
constexpr double gravity=9.80665;
int cells_number=0;
const double machine_epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
const double relative_change_in_residual=1e-6;

std::vector<std::shared_ptr<Fluid>> characterized_fluids =
    std::vector<std::shared_ptr<Fluid>>();
std::vector<std::unique_ptr<Equilibrium_Relation>> added_equilibrium_relations =
    std::vector<std::unique_ptr<Equilibrium_Relation>>();

std::shared_ptr<Mesh> mymesh;
std::shared_ptr<Rock> myrock;

void updateVariables(std::vector<std::shared_ptr<Fluid>>& _characterized_fluids){
    for(auto fluid : _characterized_fluids){
	fluid->updateProperties(term);
    };
};

double calculateAccumulation(const int _term, Fluid& _fluid, Cell& _cell, Rock& _rock){

    double past_contribution=0;
    double current_contribution=0;

    const auto _cell_index = _cell.index();
    
    for(auto equilibrium_relation : added_equilibrium_relations){
        if(equilibrium_relation->receiverFluid()->index() == _fluid.index()){
            
            const auto contributor = equilibrium_relation->contributorFluid();
            
            past_contribution = past_contribution +
                equilibrium_relation->partitionCoefficient(_term-1,_cell_index) *
                (_rock.porosity(_term-1,_cell_index) * contributor->saturation(_term-1,_cell_index)
                 / contributor->volumetricFactor(_term-1,_cell_index));

            current_contribution = current_contribution +
                equilibrium_relation->partitionCoefficient(_term,_cell_index) *
                (_rock.porosity(_term,_cell_index) * contributor->saturation(_term,_cell_index)
                 / contributor->volumetricFactor(_term,_cell_index));
            
        };
    };

    double accumulation= (_cell.getVolume()/timedelta) *
        (((_rock.porosity(_term,_cell_index)*_fluid.saturation(_term,_cell_index)
           /_fluid.volumetricFactor(_term,_cell_index)) + current_contribution)-
        ((_rock.porosity(_term-1,_cell_index)*_fluid.saturation(_term-1,_cell_index)
          /_fluid.volumetricFactor(_term-1,_cell_index)) + past_contribution));
    
    return accumulation;
}

//Change Event Name
void FluidPressureVaries(std::string& _timestamp){
    if(_timestamp == "stop"){
        updateVariables(characterized_fluids);
        
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

void calculateProperties(int& _term, int& _cell_index, int& _fluid_index){
    
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
    case 2:
        break;
    default:
        break;
    }
}

void launchPetrophysicalEngineer(){
    int option;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Rock";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        myrock = std::make_shared<Rock>(Rock());
        myrock->characterize(cells_number);
        break;
    default:
        break;
    }
};

void launchFluidsEngineer(){
    int option;
    int _dimension;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Fluid" << std::endl;
    std::cout << "2. Add Equilibrium Relation";
    
    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        characterized_fluids.push_back(std::make_shared<Fluid>(Fluid()));
	(--(characterized_fluids.end()))->get()->characterize(cells_number);
	++fluids_quantity;
        break;
    case 2:
	added_equilibrium_relations.push_back(std::make_unique<Equilibrium_Relation>());
	(--(added_equilibrium_relations.end()))->get()->add(fluids_quantity,characterized_fluids);
        break;
    default:
        break;
    }

};

void launchReservoirEngineer(){

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
