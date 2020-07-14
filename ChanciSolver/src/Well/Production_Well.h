#ifndef PRODUCTION_WELL_H
#define PRODUCTION_WELL_H

#include "Well.h"
#include "Production_Perforation.h"

class Production_Well : public Well{

private:
    
    std::vector<std::vector<double>> _rate;
    std::vector<std::vector<double>> _accumulation;    

public:

    Production_Well(const int index) : Well(index){};
    
    void perforation(Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type) override{
        
        Well::perforation(mesh, characterized_phases, type);
        insertPerforations<Production_Perforation>(mesh, characterized_phases.size());

        _rate = std::vector<std::vector<double>>(2,std::vector<double>(characterized_phases.size()));
        _accumulation = std::vector<std::vector<double>>(2,std::vector<double>(characterized_phases.size()));
        
    };

    void perforationFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type) override{
        
        Well::perforationFromFile(well_reader, mesh, characterized_phases, type);
        insertPerforationsFromFile<Production_Perforation>(well_reader, mesh, characterized_phases.size());

        _rate = std::vector<std::vector<double>>(2,std::vector<double>(characterized_phases.size()));
        _accumulation = std::vector<std::vector<double>>(2,std::vector<double>(characterized_phases.size()));
        
    };

    void updateProperties(const int& term){
        
        Well::updateProperties(term);

        if(term==1){
            _rate[1]=_rate[0];
            _accumulation[1]=_accumulation[0];
        }else{
            _rate[0]=_rate[1];
            _accumulation[0]=_accumulation[1];
        }
    };

    const std::string type() const override {return typeid(Production_Well).name();};
};



#endif /* PRODUCTION_WELL_H */
