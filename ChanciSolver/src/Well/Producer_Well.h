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
        insertPerforations<Producer_Perforate>(mesh);

        _rate = std::vector<std::vector<double>>(1,std::vector<double>(characterized_fluids.size()));
        _total_accumulated = std::vector<std::vector<double>>(1,std::vector<double>(characterized_fluids.size()));
        
    };

    void updateProperties(const int term){
        
        Well::updateProperties(term);

        _rate.push_back(_rate[term-1]);
        _total_accumulated.push_back(_total_accumulated[term-1]);
    };
};



#endif /* PRODUCER_WELL_H */
