#ifndef PRODUCER_WELL_H
#define PRODUCER_WELL_H

#include "Well.h"
#include "Producer_Perforate.h"

class Producer_Well : public Well{

 private:
    
    std::vector<std::vector<double>> _rate;
    std::vector<std::vector<double>> _total_accumulated;    

 public:

 Producer_Well(const int index) : Well(index){};
    
    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);
        insertPerforations<Producer_Perforate>(mesh, characterized_fluids.size());

        _rate = std::vector<std::vector<double>>(2,std::vector<double>(characterized_fluids.size()));
        _total_accumulated = std::vector<std::vector<double>>(2,std::vector<double>(characterized_fluids.size()));
        
    };

    void perforateFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforateFromFile(well_reader, mesh, characterized_fluids, type);
        insertPerforationsFromFile<Producer_Perforate>(well_reader, mesh, characterized_fluids.size());

        _rate = std::vector<std::vector<double>>(2,std::vector<double>(characterized_fluids.size()));
        _total_accumulated = std::vector<std::vector<double>>(2,std::vector<double>(characterized_fluids.size()));
        
    };

    void updateProperties(const int& term){
        
        Well::updateProperties(term);

        if(term==1){
            _rate[1]=_rate[0];
            _total_accumulated[1]=_total_accumulated[0];
        }else{
            _rate[0]=_rate[1];
            _total_accumulated[0]=_total_accumulated[1];
        }
    };

    const std::string type() const override {return typeid(Producer_Well).name();};
};



#endif /* PRODUCER_WELL_H */
