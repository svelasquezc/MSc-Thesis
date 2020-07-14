#include <ctime>

#include "Initial_Conditions.h"
#include "Database.h"
#include "NewtonRaphson.h"
#include "FlowFunctions.h"
#include "WellFunctions.h"

#include "VTKMesh.h"

int Phase::_main_phase_counter=0;
int Phase::_count_of_phases=0;

using BlackOilNewton = NewtonRaphson<decltype(calculateProperties), decltype(calculateFlow), decltype(calculateAccumulation),decltype(calculatePerforation),decltype(calculateWellFlow), decltype(estimateWellPressure)>;

std::unique_ptr<BlackOilNewton> my_newton;// = BlackOilNewton(0,0,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

VTKMesh vtkholder;

using namespace Database;

void insertAll(const double mytime, const int& term){
    for(auto phase : characterized_phases){
        phase->insert(vtkholder, term);
    };

    for(auto& equilibrium_relationship : added_equilibrium_relationships){
        equilibrium_relationship->insert(vtkholder, term);
    };
    
    vtkholder.write(mytime);
};

void updateVariables(const int& term){
    
    for(auto phase : characterized_phases){
	phase->updateProperties(term);
    };

    for(auto& equilibrium_relationship : added_equilibrium_relationships ){
        equilibrium_relationship->updateProperties(term);            
    };

    for(auto well : perforated_wells){
        well->updateProperties(term);
    };

    characterized_rock->updateProperties(term);
};

//Change Event Name
void PhasePressureVaries(std::string& timestamp){

    double tolerance=Initial_Conditions::divergence_tolerance*2;
    
    if(timestamp == "stop"){

        updateVariables(1);

        tolerance = my_newton->iterate(Initial_Conditions::term, *defined_mesh, perforated_wells, equations, *characterized_rock);

        while(tolerance > Initial_Conditions::relative_change_in_residual){
            updateVariables(1);
            Initial_Conditions::mytime -= Initial_Conditions::timedelta;
            Initial_Conditions::timedelta /= 2.0;
            Initial_Conditions::mytime += Initial_Conditions::timedelta;
            tolerance = my_newton->iterate(Initial_Conditions::term, *defined_mesh, perforated_wells, equations, *characterized_rock);
            
        };
        
        updateVariables(0);
        
        if(tolerance < Initial_Conditions::relative_change_in_residual/2.0){
            Initial_Conditions::timedelta *= 2.0;
        };
        
        timestamp = "continue";

        double next_time = Initial_Conditions::mytime + Initial_Conditions::timedelta;

        for(auto well : perforated_wells){
            if(next_time >= well->operativeCondition()->nextChange()){
                ++Initial_Conditions::changing_wells_quantity;
                timestamp == "change";
                well->operativeStatus(2); //Pending Change
            };
        };
    };
};

void reBuildNewton(){

    int total_equations;
    int max_non_zeros;
    int max_number_of_well_non_zeros=0;
    int well_equations=0;

    for(auto well : perforated_wells){
                
        if(well->operativeCondition()->type() == "FLOW"){
            max_number_of_well_non_zeros += well->numberOfPerforations()*2+1;
            ++well_equations;
        };                    

    };
    
    total_equations = Initial_Conditions::phases_quantity*Initial_Conditions::cells_quantity + well_equations;
    max_non_zeros = Initial_Conditions::phases_quantity*Initial_Conditions::phases_quantity*Initial_Conditions::cells_quantity + max_number_of_well_non_zeros;
    my_newton.reset();
    my_newton = std::make_unique<BlackOilNewton>(total_equations,max_non_zeros,calculateProperties,calculateFlow,calculateAccumulation,calculatePerforation,calculateWellFlow, estimateWellPressure);

    
}

void timePasses(std::string& timestamp, int& term, double& mytime, double& timedelta, double& simulationtime){
    if(timestamp == "continue" && mytime<=simulationtime){
        if(mytime == 0){
            
            for(auto cell = defined_mesh->begin(); cell != defined_mesh->end(); ++cell){
                calculateProperties(0, *cell, *characterized_rock);
            };
            
            updateVariables(1);
            
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

        //Elements which are not present in the EP,
        //but are necessary for internal management of the simulation-
        
        if(Initial_Conditions::changing_wells_quantity > 0){
            reBuildNewton();
            Initial_Conditions::changing_wells_quantity = 0;
        };
        if(mytime>=print_times[Initial_Conditions::current_print_index]){
           insertAll(mytime, term);
           ++Initial_Conditions::current_print_index;
        };
        std::cout << "Simulating interval [" << mytime << " - " << mytime + timedelta << "]" << std::endl;
        //
        
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
        defined_mesh = std::make_unique<Mesh>();
        defined_mesh->define();
        Initial_Conditions::cells_quantity = defined_mesh->getCellTotal();
        vtkholder.set(*defined_mesh);
        break;
    default:
        break;
    }
}

void launchPetrophysicalEngineer(){
    int option;
    std::unique_ptr<Interphase_Interaction> added_interphase_interaction;
    
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Rock" << std::endl;;
    std::cout << "2. Add Interphase Interaction";

    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        characterized_rock = std::make_unique<Rock>();
        characterized_rock->characterize(Initial_Conditions::cells_quantity);
        break;
    case 2:
        if(Initial_Conditions::phases_quantity >= 2){
            
            added_interphase_interaction = std::make_unique<Interphase_Interaction>();
            added_interphase_interaction->add(characterized_phases);
            
            added_interphase_interactions.push_back(std::move(added_interphase_interaction));
            
        }else{
            std::cout << "It is not possible to add an Interphase interaction with only one phase characterized."
                      << std::endl;
        }
        break;
    default:
        
        break;
    }
};

void launchPhasesEngineer(){
    int option;
    int _dimension;
    std::shared_ptr<Phase> characterized_phase;
    std::unique_ptr<Equilibrium_Relationship> added_equilibrium_relationship;
    std::cout << "Select your action" << std::endl;
    std::cout << "1. Characterize Phase" << std::endl;
    std::cout << "2. Add Equilibrium Relation";
    
    Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid option"));
    
    switch(option){
    case 1:
        characterized_phase = std::make_shared<Phase>(Initial_Conditions::phases_quantity);
        characterized_phase->characterize(Initial_Conditions::cells_quantity);
        characterized_phases.push_back(characterized_phase);
        equations.push_back(characterized_phase);
        ++Initial_Conditions::phases_quantity;
        break;
    case 2:
        if(Initial_Conditions::phases_quantity >= 2){
            added_equilibrium_relationship = std::make_unique<Equilibrium_Relationship>(Initial_Conditions::equilibrium_relationships_quantity);
            added_equilibrium_relationship->add(Initial_Conditions::cells_quantity,characterized_phases);
            added_equilibrium_relationships.push_back(std::move(added_equilibrium_relationship));
            ++Initial_Conditions::equilibrium_relationships_quantity;
            break;
        }else{
            std::cout << "It is not possible to add an Equilibrium relation with only one phase characterized."
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
    if(timestamp =="")std::cout << "1. Perforation Well" << std::endl;
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
            if(Initial_Conditions::phases_quantity >= 1){
                Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                while(type != "Production" && type != "Injection"){
                    std::cout << "Only Injection or Production wells";
                    Value_Reader::myRead(std::string("Please insert the type of well "), type, std::string("Please insert a valid input"));
                };
                ++Initial_Conditions::wells_quantity;
                if(type == "Production"){
                    well = std::make_shared<Production_Well>(Production_Well(Initial_Conditions::wells_quantity));
                }else if(type == "Injection"){
                    well = std::make_shared<Injection_Well>(Injection_Well(Initial_Conditions::wells_quantity));
                }else{
                    std::cout << "Only Injection or Production wells";
                };
            
                well->perforation(*defined_mesh, characterized_phases,type);
                perforated_wells.push_back(well);
                equations.push_back(well);
            
            }else{
                std::cout << "It is not possible to perforation wells with no phase characterized."
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
    defined_mesh->appear(timestamp,stencil);
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    PhasePressureVaries(timestamp);
    if(timestamp == "change"){
        launchReservoirEngineer(timestamp);
    }
};

void launchTriggers(std::ifstream& file_reader){
    using namespace Initial_Conditions;
    defined_mesh->appear(timestamp,stencil);
    timePasses(timestamp, term, mytime, timedelta, simulationtime);
    PhasePressureVaries(timestamp);
    reEstablishOperativeConditions(file_reader, term, timestamp);
};

void launchMenu(){
    int option;
    bool run=false;
    while(!run){
        std::cout << "Select your role or -1 for running simulation" << std::endl;
        std::cout << "1. Geomodeler"<< std::endl << "2. Petrophysical Engineer" << std::endl;
        std::cout << "3. Phases Engineer" << std::endl << "4. Reservoir Engineer";
        
        Value_Reader::myRead(std::string(""), option, std::string("Please insert a valid role"));
        
        switch(option){
        case 1:
            launchGeomodeler();
            break;
        case 2:
            launchPetrophysicalEngineer();
            break;
        case 3:
            launchPhasesEngineer();
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
    
    std::shared_ptr<Phase> characterized_phase;
    std::unique_ptr<Equilibrium_Relationship> added_equilibrium_relationship;
    std::unique_ptr<Interphase_Interaction> added_interphase_interaction;
    std::shared_ptr<Well> well;
    
    if(file_reader.is_open()){
        while(file_reader>>object){
            
            std::transform(object.begin(), object.end(),object.begin(), ::toupper);
            
            if(object == "MESH"){
                
                defined_mesh = std::make_unique<Mesh>();
                defined_mesh->defineFromFile(file_reader);
                Initial_Conditions::cells_quantity = defined_mesh->getCellTotal();
                vtkholder.set(*defined_mesh);
            }else if(object == "ROCK"){
                
                characterized_rock = std::make_unique<Rock>();
                characterized_rock->characterizeFromFile(file_reader, Initial_Conditions::cells_quantity);
            
            }else if(object == "PHASES"){
                
                file_reader>>Initial_Conditions::phases_quantity;
                
                for(int phase=0; phase<Initial_Conditions::phases_quantity;++phase){
                    
                    characterized_phase = std::make_shared<Phase>(phase);
                    characterized_phase->characterizeFromFile(file_reader, Initial_Conditions::cells_quantity);
                    characterized_phases.push_back(characterized_phase);
                    equations.push_back(characterized_phase);
                    
                };
                
            }else if(object == "EQUILIBRIUM_RELATIONSHIPS"){
                
                if(Initial_Conditions::phases_quantity >= 2){
                    file_reader>>Initial_Conditions::equilibrium_relationships_quantity;
                    for(int equilibrium_relationship=0; equilibrium_relationship<Initial_Conditions::equilibrium_relationships_quantity; ++equilibrium_relationship){
                        added_equilibrium_relationship = std::make_unique<Equilibrium_Relationship>(equilibrium_relationship);
                        added_equilibrium_relationship->addFromFile(file_reader, Initial_Conditions::cells_quantity, characterized_phases);
                        added_equilibrium_relationships.push_back(std::move(added_equilibrium_relationship));
                    };
                }else{
                    std::cout << "There is only one phase, EQUILIBRIUM_RELATIONSHIPS will be omitted\n";
                };
                
            }else if(object == "INTERPHASE_INTERACTIONS"){
                
                if(Initial_Conditions::phases_quantity >= 2){
                    file_reader>>Initial_Conditions::interphase_interactions_quantity;
                    for(int interphase_interaction=0; interphase_interaction<Initial_Conditions::interphase_interactions_quantity; ++interphase_interaction){
                        added_interphase_interaction = std::make_unique<Interphase_Interaction>(interphase_interaction);
                        added_interphase_interaction->addFromFile(file_reader, characterized_phases);
                        added_interphase_interactions.push_back(std::move(added_interphase_interaction));
                    };
                }else{
                    std::cout << "There is only one phase, INTERPHASE_INTERACTIONS will be omitted\n";
                };
                
            }else if(object == "TIME_DELTA"){
                
                file_reader>>Initial_Conditions::timedelta;
                
            }else if(object == "SIMULATION_TIME"){
                
                file_reader>>Initial_Conditions::simulationtime;

            }else if(object == "PRINT_TIMES"){

                file_reader>>Initial_Conditions::number_of_print_times;
                print_times.resize(Initial_Conditions::number_of_print_times);
                for(int print_time=0; print_time<Initial_Conditions::number_of_print_times; ++print_time){
                    file_reader>>print_times[print_time];
                };
                
            }else if(object == "WELLS"){
                if(Initial_Conditions::phases_quantity >= 1){
                    file_reader>>Initial_Conditions::wells_quantity;
                    for(int well_index=0; well_index<Initial_Conditions::wells_quantity;++well_index){
                        file_reader>>type;
                        std::transform(type.begin(), type.end(), type.begin(), ::toupper);
                        if(type == "PRODUCTION"){
                            well = std::make_shared<Production_Well>(Production_Well(well_index));
                        }else if(type == "INJECTION"){
                            well = std::make_shared<Injection_Well>(Injection_Well(well_index));
                        }else{
                            
                            std::ostringstream well_err = std::ostringstream();
                            
                            well_err<<"The type selected for well "<< well_index+1 <<": ("<< type <<") is invalid, only Production or Injection Wells"<<std::endl;
                            
                            throw std::invalid_argument(well_err.str());
                        };
                        well->perforationFromFile(file_reader,*defined_mesh, characterized_phases,type);
                        perforated_wells.push_back(well);
                        equations.push_back(well);

                    };
                }else{
                    throw std::domain_error("There are no phases characterized or WELLS Keyword appears before PHASES");
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
