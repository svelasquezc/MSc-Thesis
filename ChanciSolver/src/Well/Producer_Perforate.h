#ifndef PRODUCER_PERFORATE_H
#define PRODUCER_PERFORATE_H

#include "Perforate.h"

class Producer_Perforate : Perforate{

 private:
    
    std::vector<double> _flow;

 public:
    
    void flow(const int fluid_index, const double flow){_flow[fluid_index] = flow;};

    const double& flow(const int fluid_index) const {return _flow[fluid_index];};
};



#endif /* ṔRODUCER_PERFORATE_H */
