#ifndef EQUILIBRIUM_RELATIONSHIP_H
#define EQUILIBRIUM_RELATIONSHIP_H

#include<exception>

#include "Phase.h"

class Equilibrium_Relationship : Value_Reader{
    
private:
    int _index;
    std::weak_ptr<Phase> _contributing_phase;
    std::weak_ptr<Phase> _receiving_phase;
    std::vector<std::vector<double>>   _partition_coefficient;
    std::unique_ptr<Property_Measure> _partition_coefficient_measure;
public:
    Equilibrium_Relationship(const int index) : _index(index) {};
    Equilibrium_Relationship(){};
    void add(const int cells_quantity, std::vector<std::shared_ptr<Phase>>& characterized_phases);
    void addFromFile(std::ifstream& equilibrium_reader, const int cells_quantity, std::vector<std::shared_ptr<Phase>>& characterized_phases);
    const std::shared_ptr<Phase> contributingPhase() const {return _contributing_phase.lock();}
    const std::shared_ptr<Phase> receivingPhase()    const {return _receiving_phase.lock();}

    const double& partitionCoefficient (const int& term, const int& cell_index) const {return _partition_coefficient[term][cell_index];};             

    const std::vector<double>& partitionCoefficient (const int& term) const {return _partition_coefficient[term];};
    
    void partitionCoefficient(const int& term, const int& cell_index, const double pressure){
        _partition_coefficient[term][cell_index] =
            _partition_coefficient_measure->interpolate(pressure);
    };
    
    void updateProperties(const int& term){
        if(term==1){
            _partition_coefficient[1]=_partition_coefficient[0];
        }else{
            _partition_coefficient[0]=_partition_coefficient[1];
        }
    };

    template<typename VTKType> void insert(VTKType& vtkcontext, const int term){
        
        std::ostringstream fullname;

        fullname << receivingPhase()->type() << " in " << contributingPhase()->type() << " Ratio";
        
        vtkcontext.appendScalar(fullname.str(), _partition_coefficient[term]);
    };
};

void Equilibrium_Relationship::add(const int cells_quantity, std::vector<std::shared_ptr<Phase>>& characterized_phases){

    int counter;
    int contributing;
    int receiving;

    std::stringstream ref_value = std::stringstream();
    
    if(characterized_phases.size() >= 2){

        _partition_coefficient = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
        
        std::cout << "Please select the contributing Phase: " << std::endl;
        for (counter=0; counter<characterized_phases.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
        };

        while(true){
            myRead(std::string(""), contributing, std::string("Please insert a valid index"));
            if(contributing>0 && contributing<=characterized_phases.size()){
                _contributing_phase = characterized_phases[contributing-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };
        
        std::cout << "Select the receiving Phase: " << std::endl;
        counter = 0;
        for (counter=0; counter<characterized_phases.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
        };
        while(true){
            myRead(std::string(""), receiving, std::string("Please insert a valid index"));
            if(receiving>0 && receiving<=characterized_phases.size()){
                _receiving_phase = characterized_phases[receiving-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
        };

        ref_value << _receiving_phase.lock()->print() << " in " << _contributing_phase.lock()->print() << " ratio";
        
        _partition_coefficient_measure =
            std::make_unique<Property_Measure>(Property_Measure(std::string("Pressure"),ref_value.str()));
        
        _partition_coefficient_measure->readMe();
        
    }else{
        
        std::cout << "It is not possible to add an Equilibrium relation with only one phase characterized."
                  << std::endl;
        
    };
    
};

void Equilibrium_Relationship::addFromFile(std::ifstream& equilibrium_reader, const int cells_quantity, std::vector<std::shared_ptr<Phase>>& characterized_phases){
    
    std::string element;
    int contributing;
    int receiving;
    std::stringstream ref_value = std::stringstream();
        
    if(characterized_phases.size() >= 2){
        
        _partition_coefficient = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
        
        while(equilibrium_reader >> element){

            std::transform(element.begin(), element.end(),element.begin(), ::toupper);

            if(element == "CONTRIBUTING_PHASE"){
                
                equilibrium_reader >> contributing;

                if(contributing < 1 && contributing>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Contributing phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                }else{
                    _contributing_phase = characterized_phases[contributing-1];
                };
                    
            }else if(element == "RECEIVING_PHASE"){
                
                equilibrium_reader >> receiving;

                if(receiving < 1 && receiving>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Receiving phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                }else{
                    _receiving_phase = characterized_phases[receiving-1];
                };
                
            }else if(element == "PARTITION_COEFFICIENT"){

                ref_value << _receiving_phase.lock()->print() << " in " << _contributing_phase.lock()->print() << " ratio";
        
                _partition_coefficient_measure =
                    std::make_unique<Property_Measure>(Property_Measure(std::string("Pressure"),ref_value.str()));
                _partition_coefficient_measure->readFromFile(equilibrium_reader);
                break;
            };
            
        };
        
    };
};

#endif /* EQUILIBRIUM_RELATIONSHIP_H */
