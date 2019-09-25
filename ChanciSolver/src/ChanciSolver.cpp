#include <ctime>

#include "Initial_Conditions.h"
#include "Database.h"
#include "NewtonRaphson.h"
#include "FlowFunctions.h"
#include "WellFunctions.h"

#include "VTKMesh.h"

int Fluid::_count_of_principals=0;
int Fluid::_count_of_fluids=0;

using BlackOilNewton = NewtonRaphson<decltype(calculateProperties), decltype(calculateFlow), decltype(calculateAccumulation),decltype(calculatePerforation),decltype(calculateWellFlow), decltype(estimateWellPressure)>;

std::unique_ptr<BlackOilNewton> my_newton;// = BlackOilNewton(0,0,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

VTKMesh vtkholder;

using namespace Database;

void insertAll(const double mytime, const int& term){
    for(auto fluid : characterized_fluids){
        fluid->insert(vtkholder, term);
    };

    for(auto& equilibrium_relation : added_equilibrium_relations){
        equilibrium_relation->insert(vtkholder, term);
    };
    
    vtkholder.write(mytime);
};

void updateVariables(const int& term, Rock& rock){
    
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

//Change Event Name
void FluidPressureVaries(std::string& timestamp){
    if(timestamp == "stop"){

        updateVariables(1, *myrock);

        my_newton->iterate(Initial_Conditions::term, *mymesh, perforated_wells, equations, *myrock);

        updateVariables(0, *myrock);
        
        timestamp = "continue";

        double next_time = Initial_Conditions::mytime + Initial_Conditions::timedelta;

        for(auto well : perforated_wells){
            if(next_time >= well->operativeCondition()->nextChange()){
                ++Initial_Conditions::changing_wells;
                timestamp == "change";
                well->operativeStatus(2); //Pending Change
            };
        };
    }
};

void reBuildNewton(){

    int total_equations;
    int max_non_zeros;
    int max_number_of_well_non_zeros=0;
    int well_equations=0;

    for(auto well : perforated_wells){
                
        if(well->operativeCondition()->type() == "FLOW"){
            max_number_of_well_non_zeros += well->numberOfPerforates()*2+1;
            ++well_equations;
        };                    

    };
    
    total_equations = Initial_Conditions::fluids_quantity*Initial_Conditions::cells_number + well_equations;
    max_non_zeros = Initial_Conditions::fluids_quantity*Initial_Conditions::fluids_quantity*Initial_Conditions::cells_number + max_number_of_well_non_zeros;
    my_newton.reset();
    my_newton = std::make_unique<BlackOilNewton>(total_equations,max_non_zeros,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

    
}

void timePasses(std::string& timestamp, int& term, double& mytime, double& timedelta, double& simulationtime){
    if(timestamp == "continue" && mytime<=simulationtime){
        if(mytime == 0){
            for(auto cell = mymesh->begin(); cell != mymesh->end(); ++cell){
                calculateProperties(0, *cell, *myrock);
            };

            for(auto well : perforated_wells){
                
                if(well->operativeStatus() == 1){ //Just Changed

                    if(well->operativeCondition()->type() == "FLOW"){
                        estimateWellPressure(0, well);
                    }else{
                        calculateWellFlow(0, well);
                    };
                    
                    well->operativeStatus(0); // Stable
                };
            };
            
            reBuildNewton();
            
            ++term;
            
        };

        if(Initial_Conditions::changing_wells > 0){
            reBuildNewton();
            Initial_Conditions::changing_wells = 0;
        };

        insertAll(mytime, term);
        
        std::cout << "Simulating interval [" << mytime << " - " << mytime + timedelta << "]" << std::endl;
        
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
        Initial_Conditions::cells_number = mymesh->getCellTotal();
        vtkholder.set(*mymesh);
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
        myrock->characterize(Initial_Conditions::cells_number);
        break;
    case 2:
        if(Initial_Conditions::fluids_quantity >= 2){
            
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
        characterized_fluid = std::make_shared<Fluid>(Initial_Conditions::fluids_quantity);
        characterized_fluid->characterize(Initial_Conditions::cells_number);
        characterized_fluids.push_back(characterized_fluid);
        equations.push_back(characterized_fluid);
        ++Initial_Conditions::fluids_quantity;
        break;
    case 2:
        if(Initial_Conditions::fluids_quantity >= 2){
            added_equilibrium_relation = std::make_unique<Equilibrium_Relation>(Initial_Conditions::equilibrium_relations_quantity);
            added_equilibrium_relation->add(Initial_Conditions::cells_number,characterized_fluids);
            added_equilibrium_relations.push_back(std::move(added_equilibrium_relation));
            ++Initial_Conditions::equilibrium_relations_quantity;
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
            if(Initial_Conditions::fluids_quantity >= 1){
                Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                while(type != "Producer" && type != "Injector"){
                    std::cout << "Only Injector or Producer wells";
                    Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                };
                ++Initial_Conditions::wells_quantity;
                if(type == "Producer"){
                    well = std::make_shared<Producer_Well>(Producer_Well(Initial_Conditions::wells_quantity));
                }else if(type == "Injector"){
                    well = std::make_shared<Injector_Well>(Injector_Well(Initial_Conditions::wells_quantity));
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
        while(index<1 && index>=Initial_Conditions::wells_quantity){
            Value_Reader::myRead(std::string("Please select an index between the range "), index, std::string("Please insert a valid input"));
        };

        perforated_wells[index-1]->establish(Initial_Conditions::term, Initial_Conditions::timestamp);
        
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
    using namespace Initial_Conditions;
    mymesh->appear(timestamp,stencil);
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    FluidPressureVaries(timestamp);
    if(timestamp == "change"){
        launchReservoirEngineer(timestamp);
    }
};

void launchTriggers(std::ifstream& file_reader){
    using namespace Initial_Conditions;
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
            launchReservoirEngineer(Initial_Conditions::timestamp);
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
                Initial_Conditions::cells_number = mymesh->getCellTotal();
                vtkholder.set(*mymesh);
            }else if(object == "ROCK"){
                
                myrock = std::make_unique<Rock>();
                myrock->characterizeFromFile(file_reader, Initial_Conditions::cells_number);
            
            }else if(object == "FLUIDS"){
                
                file_reader>>Initial_Conditions::fluids_quantity;
                
                for(int fluid=0; fluid<Initial_Conditions::fluids_quantity;++fluid){
                    
                    characterized_fluid = std::make_shared<Fluid>(fluid);
                    characterized_fluid->characterizeFromFile(file_reader, Initial_Conditions::cells_number);
                    characterized_fluids.push_back(characterized_fluid);
                    equations.push_back(characterized_fluid);
                    
                };
                
            }else if(object == "EQUILIBRIUM_RELATIONS"){
                
                if(Initial_Conditions::fluids_quantity >= 2){
                    file_reader>>Initial_Conditions::equilibrium_relations_quantity;
                    for(int equilibrium_relation=0; equilibrium_relation<Initial_Conditions::equilibrium_relations_quantity; ++equilibrium_relation){
                        added_equilibrium_relation = std::make_unique<Equilibrium_Relation>(equilibrium_relation);
                        added_equilibrium_relation->addFromFile(file_reader, Initial_Conditions::cells_number, characterized_fluids);
                        added_equilibrium_relations.push_back(std::move(added_equilibrium_relation));
                    };
                }else{
                    std::cout << "There is only one fluid, EQUILIBRIUM_RELATIONS will be omitted\n";
                };
                
            }else if(object == "INTERFLUID_INTERACTIONS"){
                
                if(Initial_Conditions::fluids_quantity >= 2){
                    file_reader>>Initial_Conditions::interfluid_interactions_quantity;
                    for(int interfluid_interaction=0; interfluid_interaction<Initial_Conditions::interfluid_interactions_quantity; ++interfluid_interaction){
                        added_interfluid_interaction = std::make_unique<Interfluid_Interaction>(interfluid_interaction);
                        added_interfluid_interaction->addFromFile(file_reader, characterized_fluids);
                        added_interfluid_interactions.push_back(std::move(added_interfluid_interaction));
                    };
                }else{
                    std::cout << "There is only one fluid, INTERFLUID_INTERACTIONS will be omitted\n";
                };
                
            }else if(object == "TIME_DELTA"){
                
                file_reader>>Initial_Conditions::timedelta;
                
            }else if(object == "SIMULATION_TIME"){
                
                file_reader>>Initial_Conditions::simulationtime;
                
            }else if(object == "WELLS"){
                if(Initial_Conditions::fluids_quantity >= 1){
                    file_reader>>Initial_Conditions::wells_quantity;
                    for(int well_index=0; well_index<Initial_Conditions::wells_quantity;++well_index){
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
                if(Initial_Conditions::wells_quantity >= 1){
                    for(auto well : perforated_wells){
                        file_reader >> object;
                        std::transform(object.begin(), object.end(),object.begin(), ::toupper);
                        if(object == "WELL"){
                            int index;
                            file_reader >> index;
                            if(well->index() == index - 1){
                                well->establishFromFile(file_reader, Initial_Conditions::term, Initial_Conditions::timestamp);
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
        
        Initial_Conditions::timestamp="continue";
        tstart = std::time(0);
        while(Initial_Conditions::mytime<Initial_Conditions::simulationtime){
            launchTriggers();
        };
        tend = std::time(0);
        
        std::cout << "Elapsed time: "<< std::difftime(tend, tstart) <<" seconds(s)."<<std::endl;
        
    }else{
        std::ifstream file_reader(argv[1], std::ifstream::in);
        launchFromFile(file_reader);
        Initial_Conditions::timestamp="continue";
        tstart = std::time(0);
        while(Initial_Conditions::mytime<Initial_Conditions::simulationtime){
            launchTriggers(file_reader);
        };
        tend = std::time(0);

        std::cout << "Elapsed time: "<< std::difftime(tend, tstart) <<" seconds(s)."<<std::endl;
    };

    insertAll(Initial_Conditions::mytime, Initial_Conditions::term);
    
    return 0;
};
