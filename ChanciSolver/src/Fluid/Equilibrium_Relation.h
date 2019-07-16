#ifndef EQUILIBRIUM_RELATION_H
#define EQUILIBRIUM_RELATION_H

#include "Fluid.h"

class Equilibrium_Relation : Value_Reader{
    
 private:
    int index;
    std::shared_ptr<Fluid> _contributor_fluid;
    std::shared_ptr<Fluid> _receiver_fluid;
    std::vector<std::vector<double>>   _partition_coefficient;
    std::unique_ptr<Measured_Property> _measured_partition_coefficient;
 public:
    Equilibrium_Relation(){};
    void add(const int fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids);
    const std::shared_ptr<Fluid>& contributorFluid() const {return _contributor_fluid;}
    const std::shared_ptr<Fluid>& receiverFluid()    const {return _receiver_fluid;}

    const double& partitionCoefficient (const int& term, const int& cell_index) const {return _partition_coefficient[term][cell_index];};             

    void partitionCoefficient(int& term, int& cell_index){
        _partition_coefficient[term][cell_index] =
            _measured_partition_coefficient->interpolate(_contributor_fluid->pressure(term,cell_index));
    };
    
    void updateProperties(const int& term){
        _partition_coefficient.push_back(_partition_coefficient[term-1]);
    }
};

void Equilibrium_Relation::add(const int fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids){

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
            if(contributor>0 && contributor<=MyFluids.size()){
                _contributor_fluid = MyFluids[contributor-1];
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
                _receiver_fluid = MyFluids[receiver-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
        };

        ref_value << _receiver_fluid->print() << " in " << _contributor_fluid->print() << " ratio";
        
        _measured_partition_coefficient =
            std::make_unique<Measured_Property>(Measured_Property(std::string("Pressure"),ref_value.str()));
        
        _measured_partition_coefficient->readMe();
        
    }else{
        
        std::cout << "It is not possible to add a Equilibrium relation with only one fluid characterized."
                  << std::endl;
        
    };
    
};


#endif /* EQUILIBRIUM_RELATION_H */
