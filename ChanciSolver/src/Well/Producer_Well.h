#ifndef PRODUCER_WELL_H
#define PRODUCER_WELL_H

#include "Well.h"
#include "Producer_Perforate.h"

class Producer_Well : public Well, public std::enable_shared_from_this<Producer_Well>{

 private:
    
    std::vector<std::vector<double>> _rate;
    std::vector<std::vector<double>> _total_accumulated;    

 public:

 Producer_Well() : Well(){};
    
    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);
        insertPerforations<Producer_Perforate>(mesh);

        _rate = std::vector<std::vector<double>>();
        _total_accumulated = std::vector<std::vector<double>>();
        
    };
};



#endif /* PRODUCER_WELL_H */
