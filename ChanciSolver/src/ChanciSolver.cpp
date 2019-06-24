
#include "Mesh.h"

std::string timestamp="";
double simulationtime=0;
int fluids_quantity=0;
int stencil[2] = {-1,1};
int equilibrium_relations_quantity=0;
constexpr double gravity=9.80665;
int cells_quantity=0;
const double machine_epsilon = 1e-16;
const double relative_change_in_solution=1e-6;
Mesh mymesh;
//Rock
void launchTriggers(){
  mymesh.appear(timestamp,stencil);
};
void launchGeomodeler(){
    int option;
    int _dimension;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Define Mesh" << std::endl;
    std::cin >> option;
    switch(option){
    case 1:
        mymesh = Mesh();
        mymesh.defineMesh();
        break;
    case 2:
        break;
    default:
        break;
    }
    
}

void launchMenu(){
    int option;
    while(true){
        std::cout << "Select your role" << std::endl;
        std::cout << "1. Geomodeler"<< std::endl << "2. Petrophysical Engineer" << std::endl;
        std::cout << "3. Fluids Engineer" << std::endl << "4. Reservoir Engineer" << std::endl;
        try{
            std::cin >> option;
            switch(option){
            case 1:
                launchGeomodeler();
		break;
            }
	    launchTriggers();
	    break;
        }catch(std::exception e){
            
            continue;
        }
    }
}

int main(){
    launchMenu();
    timestamp="s";
    launchTriggers();
    return 0;
};
