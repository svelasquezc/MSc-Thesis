
#include "Mesh.h"
#include "Rock.h"
#include "Equilibrium_Relation.h"
#include "Jacobian.h"

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
const double machine_epsilon = 1e-16;
const double relative_change_in_solution=1e-6;

std::vector<std::shared_ptr<Fluid>> characterized_fluids =
    std::vector<std::shared_ptr<Fluid>>();
std::vector<std::unique_ptr<Equilibrium_Relation>> added_equilibrium_relations =
    std::vector<std::unique_ptr<Equilibrium_Relation>>();

Mesh mymesh;
Rock myrock;




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
    mymesh.appear(timestamp,stencil);
    
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
        mymesh = Mesh();
        mymesh.defineMesh();
        cells_number = mymesh.getCellTotal();
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
        myrock = Rock();
        myrock.characterize(cells_number);
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
