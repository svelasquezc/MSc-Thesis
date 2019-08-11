#ifndef INJECTOR_WELL_H
#define INJECTOR_WELL_H

#include "Fluid.h"

class Injector_Well : Well{

 private:

    std::shared_ptr<Fluid> _injection_fluid;
    
    std::vector<std::vector<double>> _rate;
    std::vector<std::vector<double>> _total_accumulated;    

 public:
    
};



#endif /* INJECTOR_WELL_H */
h