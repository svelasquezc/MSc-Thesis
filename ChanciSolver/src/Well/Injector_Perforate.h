#ifndef INJECTOR_PERFORATE_H
#define INJECTOR_PERFORATE_H

#include "Perforate.h"

class Injector_Perforate : public Perforate{

 private:
    
    double _flow=0;

 public:

 Injector_Perforate() : Perforate(){};

    void flow(double flow){_flow=flow;};
    
    const double totalFlow() const override {return _flow;};
    const double& flow() const{return _flow;};
};



#endif /* INJECTOR_PERFORATE_H */
