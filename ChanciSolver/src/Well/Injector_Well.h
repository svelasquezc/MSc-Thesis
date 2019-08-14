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

    void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type) override{
        
        Well::perforate(mesh, characterized_fluids, type);

        insertPerforations<Injector_Perforate>(mesh);
        
    };
};


#endif /* INJECTOR_WELL_H */
