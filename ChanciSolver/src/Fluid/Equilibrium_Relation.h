#ifndef EQUILIBRIUM_RELATION_H
#define EQUILIBRIUM_RELATION_H

#include "Fluid.h"

class Equilibrium_Relation{
    
 private:
    int index;
    std::shared_ptr<Fluid> contributor_fluid;
    std::shared_ptr<Fluid> receiver_fluid;
    std::vector<std::vector<double>> partition_coefficient;
    std::unique_ptr<Measured_Property> measured_partition_coefficient;
 public:
    Equilibrium_Relation(){
        measured_partition_coefficient = std::make_unique<Measured_Property>();
    };
    void add(const int& fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids);
    
};

void Equilibrium_Relation::add(const int& fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids){

    int counter;
    int contributor;
    int receiver;
    
    if(fluids_quantity >= 2){
        
        
        std::cout << "Please select the contributor Fluid: " << std::endl;
        for (counter=0; counter<MyFluids.size(); ++counter){
            std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
        };

        while(true){
            try{
                std::cin >> contributor;
                if(contributor>0 && contributor<=MyFluids.size()){
                    contributor_fluid = MyFluids[contributor-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            }catch(std::exception e){
                std::cout << "Please insert a valid index" << std::endl;
            }
        };
        
        std::cout << "Select the receiver Fluid: " << std::endl;
        counter = 0;
        for (counter=0; counter<MyFluids.size(); ++counter){
            std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
        };
        while(true){
            try{
                std::cin >> receiver;
                if(receiver>0 && receiver<=MyFluids.size()){
                    receiver_fluid = MyFluids[receiver-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            }catch(std::exception e){
                std::cout << "Please insert a valid index" << std::endl;
            }
        };

        measured_partition_coefficient->readMe();
        
    }else{
        
        std::cout << "It is not possible to add a Equilibrium relation with only one fluid characterized."
                  << std::endl;
        
    };
    
}

#endif /* EQUILIBRIUM_RELATION_H */
