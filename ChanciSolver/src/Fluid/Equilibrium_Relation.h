#ifndef EQUILIBRIUM_RELATION_H
#define EQUILIBRIUM_RELATION_H

#include<exception>

#include "Fluid.h"

class Equilibrium_Relation : Value_Reader{
    
 private:
    int _index;
    std::weak_ptr<Fluid> _contributor_fluid;
    std::weak_ptr<Fluid> _receiver_fluid;
    std::vector<std::vector<double>>   _partition_coefficient;
    std::unique_ptr<Measured_Property> _measured_partition_coefficient;
 public:
 Equilibrium_Relation(const int index) : _index(index) {};
    Equilibrium_Relation(){};
    void add(const int cells_number, std::vector<std::shared_ptr<Fluid>>& MyFluids);
    void addFromFile(std::ifstream& equilibrium_reader, const int cells_number, std::vector<std::shared_ptr<Fluid>>& MyFluids);
    const std::shared_ptr<Fluid> contributorFluid() const {return _contributor_fluid.lock();}
    const std::shared_ptr<Fluid> receiverFluid()    const {return _receiver_fluid.lock();}

    const double& partitionCoefficient (const int& term, const int& cell_index) const {return _partition_coefficient[term][cell_index];};             

    void partitionCoefficient(const int& term, const int& cell_index, const double pressure){
        _partition_coefficient[term][cell_index] =
            _measured_partition_coefficient->interpolate(pressure);
    };
    
    void updateProperties(const int& term){
        _partition_coefficient.push_back(_partition_coefficient[term-1]);
    }
};

void Equilibrium_Relation::add(const int cells_number, std::vector<std::shared_ptr<Fluid>>& characterized_fluids){

    int counter;
    int contributor;
    int receiver;

    std::stringstream ref_value = std::stringstream();
    
    if(characterized_fluids.size() >= 2){

        _partition_coefficient = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
        
        std::cout << "Please select the contributor Fluid: " << std::endl;
        for (counter=0; counter<characterized_fluids.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
        };

        while(true){
            myRead(std::string(""), contributor, std::string("Please insert a valid index"));
            if(contributor>0 && contributor<=characterized_fluids.size()){
                _contributor_fluid = characterized_fluids[contributor-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };
        
        std::cout << "Select the receiver Fluid: " << std::endl;
        counter = 0;
        for (counter=0; counter<characterized_fluids.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
        };
        while(true){
            myRead(std::string(""), receiver, std::string("Please insert a valid index"));
            if(receiver>0 && receiver<=characterized_fluids.size()){
                _receiver_fluid = characterized_fluids[receiver-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
        };

        ref_value << _receiver_fluid.lock()->print() << " in " << _contributor_fluid.lock()->print() << " ratio";
        
        _measured_partition_coefficient =
            std::make_unique<Measured_Property>(Measured_Property(std::string("Pressure"),ref_value.str()));
        
        _measured_partition_coefficient->readMe();
        
    }else{
        
        std::cout << "It is not possible to add an Equilibrium relation with only one fluid characterized."
                  << std::endl;
        
    };
    
};

void Equilibrium_Relation::addFromFile(std::ifstream& equilibrium_reader, const int cells_number, std::vector<std::shared_ptr<Fluid>>& characterized_fluids){
    
    std::string element;
    int contributor;
    int receiver;
    std::stringstream ref_value = std::stringstream();
        
    if(characterized_fluids.size() >= 2){
        
        _partition_coefficient = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    
        while(equilibrium_reader >> element){

            std::transform(element.begin(), element.end(),element.begin(), ::toupper);

            if(element == "CONTRIBUTOR_FLUID"){
                
                equilibrium_reader >> contributor;

                if(contributor < 1 && contributor>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Contributor fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                }else{
                    _contributor_fluid = characterized_fluids[contributor-1];
                };
                    
            }else if(element == "RECEIVER_FLUID"){
                
                equilibrium_reader >> receiver;

                if(receiver < 1 && receiver>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Receiver fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                }else{
                    _receiver_fluid = characterized_fluids[receiver-1];
                };
                
            }else if(element == "PARTITION_COEFFICIENT"){

                ref_value << _receiver_fluid.lock()->print() << " in " << _contributor_fluid.lock()->print() << " ratio";
        
                _measured_partition_coefficient =
                    std::make_unique<Measured_Property>(Measured_Property(std::string("Pressure"),ref_value.str()));
                _measured_partition_coefficient->readFromFile(equilibrium_reader);
                break;
            };
            
        };
        
    };
};

#endif /* EQUILIBRIUM_RELATION_H */
