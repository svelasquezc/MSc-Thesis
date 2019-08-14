#ifndef INJECTOR_PERFORATE_H
#define INJECTOR_PERFORATE_H

#include "Perforate.h"

class Injector_Perforate : public Perforate{

 private:
    
    double _flow;

 public:

    void flow(const double flow) {_flow = flow;};
    const double& flow() const{return _flow;};
};



#endif /* INJECTOR_PERFORATE_H */
