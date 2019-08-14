#ifndef INJECTOR_WELL_H
#define INJECTOR_WELL_H

#include "Well.h"
#include "Fluid.h"


class Injector_Well : public Well{

 private:

    std::shared_ptr<Fluid> _injection_fluid;
    
    std::vector<double> _rate;
    std::vector<double> _total_accumulated;    

 public:
    
};


#endif /* INJECTOR_WELL_H */
