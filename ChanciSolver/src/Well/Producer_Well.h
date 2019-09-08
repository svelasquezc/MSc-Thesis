#ifndef PRODUCER_WELL_H
#define PRODUCER_WELL_H

#include "Well.h"
#include "Producer_Perforate.h"

class Producer_Well : public Well{

 private:
    
    std::map<const std::string, std::vector<double>> _rate;
    std::map<const std::string, std::vector<double>> _total_accumulated;    

 public:

 Producer_Well(const int index) : Well(index){};
    
    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);
        insertPerforations<Producer_Perforate>(mesh, characterized_fluids.size());

        _rate["N"] = std::vector<double>(characterized_fluids.size());
        _total_accumulated["N"] = std::vector<double>(characterized_fluids.size());

        _rate["K"] = std::vector<double>(characterized_fluids.size());
        _total_accumulated["K"] = std::vector<double>(characterized_fluids.size());
        
    };

    void perforateFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforateFromFile(well_reader, mesh, characterized_fluids, type);
        insertPerforationsFromFile<Producer_Perforate>(well_reader, mesh, characterized_fluids.size());

        _rate["N"] = std::vector<double>(characterized_fluids.size());
        _total_accumulated["N"] = std::vector<double>(characterized_fluids.size());

        _rate["K"] = std::vector<double>(characterized_fluids.size());
        _total_accumulated["K"] = std::vector<double>(characterized_fluids.size());
        
    };

    void updateProperties(const std::string& term){
        
        Well::updateProperties(term);

        if(term=="K"){
            _rate["K"]=_rate["N"];
            _total_accumulated["K"]=_total_accumulated["N"];
        }else{
            _rate["N"]=_rate["K"];
            _total_accumulated["N"]=_total_accumulated["K"];
        }
    };

    const std::string type() const override {return typeid(Producer_Well).name();};
};



#endif /* PRODUCER_WELL_H */
