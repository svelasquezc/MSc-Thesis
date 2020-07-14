#ifndef INJECTION_PERFORATION_H
#define INJECTION_PERFORATION_H

#include "Perforation.h"

class Injection_Perforation : public Perforation{

 private:
    
    double _flow=0;

 public:

 Injection_Perforation() : Perforation(){};
    
 Injection_Perforation(const int& fluids_quantity) : Perforation(){};

    void flow(double flow){_flow=flow;};
    
    const double totalFlow() const override {return _flow;};
    const double& flow() const{return _flow;};

    const std::string type() const override {return typeid(Injection_Perforation).name();};
};



#endif /* INJECTION_PERFORATION_H */
