#ifndef INJECTION_WELL_H
#define INJECTION_WELL_H

#include "Well.h"
#include "Phase.h"
#include "Injection_Perforation.h"

class Injection_Well : public Well{

 private:

    std::weak_ptr<Phase> _injection_phase;
    
    std::vector<double> _rate;
    std::vector<double> _accumulation;    

 public:

 Injection_Well(const int index) : Well(index){};
    
    void perforation(Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type) override{
        
        Well::perforation(mesh, characterized_phases, type);

        insertPerforations<Injection_Perforation>(mesh, characterized_phases.size());

        int injection;
        
        std::cout << "Please select the injection Phase: " << std::endl;
        for (int counter=0; counter<characterized_phases.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
        };

        while(true){
            Value_Reader::myRead(std::string(""), injection, std::string("Please insert a valid index"));
            if(injection>0 && injection<=characterized_phases.size()){
                _injection_phase = characterized_phases[injection-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };

        _rate = std::vector<double>(2,0.0);
        _accumulation = std::vector<double>(2,0.0);
        
    };

    void perforationFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type) override{
        
        Well::perforationFromFile(well_reader, mesh, characterized_phases, type);

        insertPerforationsFromFile<Injection_Perforation>(well_reader, mesh, characterized_phases.size());

        int injection;
        std::string element;

        while(well_reader >> element){
            
            std::transform(element.begin(), element.end(),element.begin(), ::toupper);
            
            if(element == "INJECTION_PHASE"){
                well_reader >> injection;

                if(injection < 1 && injection>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Injection phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                }else{
                    _injection_phase = characterized_phases[injection-1];
                };
                
                break;
            };
        };

        _rate = std::vector<double>(2,0.0);
        _accumulation = std::vector<double>(2,0.0);
        
    };

    const std::shared_ptr<Phase> injectionPhase() const {return  _injection_phase.lock();};

    void updateProperties(const int& term){
        Well::updateProperties(term);

        if(term==1){
            _rate[1]=_rate[0];
            _accumulation[1]=_accumulation[0];
        }else{
            _rate[1]=_rate[0];
            _accumulation[1]=_accumulation[0];
        }
    };

    const std::string type() const override {return typeid(Injection_Well).name();};
};


#endif /* INJECTION_WELL_H */
