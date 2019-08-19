#ifndef INJECTOR_WELL_H
#define INJECTOR_WELL_H

#include "Well.h"
#include "Fluid.h"
#include "Injector_Perforate.h"

class Injector_Well : public Well{

 private:

    std::shared_ptr<Fluid> _injection_fluid;
    
    std::vector<double> _rate;
    std::vector<double> _total_accumulated;    

 public:

 Injector_Well() : Well(){};
    
    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);

        insertPerforations<Injector_Perforate>(mesh);

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

        _rate=std::vector<double>();
        _total_accumulated=std::vector<double>();
        
    };

    const std::shared_ptr<Fluid>& injectionFluid() const {return  _injection_fluid;};

    void updateProperties(const int term){
        Well::updateProperties(term);

        _rate.push_back(_rate[term-1]);
        _total_accumulated.push_back(_total_accumulated[term-1]);
    };
};


#endif /* INJECTOR_WELL_H */
