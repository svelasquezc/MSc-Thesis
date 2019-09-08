#ifndef INJECTOR_WELL_H
#define INJECTOR_WELL_H

#include "Well.h"
#include "Fluid.h"
#include "Injector_Perforate.h"

class Injector_Well : public Well{

 private:

    std::weak_ptr<Fluid> _injection_fluid;
    
    std::map<const std::string, double> _rate;
    std::map<const std::string, double> _total_accumulated;    

 public:

 Injector_Well(const int index) : Well(index){};
    
    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);

        insertPerforations<Injector_Perforate>(mesh, characterized_fluids.size());

        int injection;
        
        std::cout << "Please select the injection Fluid: " << std::endl;
        for (int counter=0; counter<characterized_fluids.size(); ++counter){
            std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
        };

        while(true){
            Value_Reader::myRead(std::string(""), injection, std::string("Please insert a valid index"));
            if(injection>0 && injection<=characterized_fluids.size()){
                _injection_fluid = characterized_fluids[injection-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };

        _rate ["N"]=0.0;
        _total_accumulated["N"]= 0.0;

        _rate ["K"]=0.0;
        _total_accumulated["K"]= 0.0;
        
    };

    void perforateFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforateFromFile(well_reader, mesh, characterized_fluids, type);

        insertPerforationsFromFile<Injector_Perforate>(well_reader, mesh, characterized_fluids.size());

        int injection;
        std::string element;

        while(well_reader >> element){
            
            std::transform(element.begin(), element.end(),element.begin(), ::toupper);
            
            if(element == "INJECTION_FLUID"){
                well_reader >> injection;

                if(injection < 1 && injection>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Injection fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                }else{
                    _injection_fluid = characterized_fluids[injection-1];
                };
                
                break;
            };
        };

        _rate ["N"]=0.0;
        _total_accumulated["N"]= 0.0;
        
        _rate ["K"]=0.0;
        _total_accumulated["K"]= 0.0;
        
    };

    const std::shared_ptr<Fluid> injectionFluid() const {return  _injection_fluid.lock();};

    void updateProperties(const std::string& term){
        Well::updateProperties(term);

        if(term=="K"){
            _rate["K"]=_rate["N"];
            _total_accumulated["K"]=_total_accumulated["N"];
        }else{
            _rate["K"]=_rate["N"];
            _total_accumulated["K"]=_total_accumulated["N"];
        }
    };

    const std::string type() const override {return typeid(Injector_Well).name();};
};


#endif /* INJECTOR_WELL_H */
