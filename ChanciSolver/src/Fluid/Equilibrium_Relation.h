#ifndef EQUILIBRIUM_RELATION_H
#define EQUILIBRIUM_RELATION_H

#include "Fluid.h"

class Equilibrium_Relation : Value_Reader{
    
 private:
    int index;
    std::shared_ptr<Fluid> contributor_fluid;
    std::shared_ptr<Fluid> receiver_fluid;
    std::vector<std::vector<double>> partition_coefficient;
    std::unique_ptr<Measured_Property> measured_partition_coefficient;
 public:
    Equilibrium_Relation(){
        
    };
    void add(const int& fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids);
    
};

void Equilibrium_Relation::add(const int& fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids){

    int counter;
    int contributor;
    int receiver;

    std::stringstream ref_value = std::stringstream();
    
    if(fluids_quantity >= 2){
        
        
        std::cout << "Please select the contributor Fluid: " << std::endl;
        for (counter=0; counter<MyFluids.size(); ++counter){
            std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
        };

        while(true){
            myRead(std::string(""), contributor, std::string("Please insert a valid index"));
            std::cin >> contributor;
            if(contributor>0 && contributor<=MyFluids.size()){
                contributor_fluid = MyFluids[contributor-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };
        
        std::cout << "Select the receiver Fluid: " << std::endl;
        counter = 0;
        for (counter=0; counter<MyFluids.size(); ++counter){
            std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
        };
        while(true){
            myRead(std::string(""), receiver, std::string("Please insert a valid index"));
            if(receiver>0 && receiver<=MyFluids.size()){
                receiver_fluid = MyFluids[receiver-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
        };

        ref_value << receiver_fluid->print() << " in " << contributor_fluid->print() << " ratio";
        
        measured_partition_coefficient =
            std::make_unique<Measured_Property>(Measured_Property(std::string("Pressure"),ref_value.str()));
        
        measured_partition_coefficient->readMe();
        
    }else{
        
        std::cout << "It is not possible to add a Equilibrium relation with only one fluid characterized."
                  << std::endl;
        
    };
    
}

#endif /* EQUILIBRIUM_RELATION_H */
