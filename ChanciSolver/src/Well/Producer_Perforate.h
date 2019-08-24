#ifndef PRODUCER_PERFORATE_H
#define PRODUCER_PERFORATE_H

#include "Perforate.h"
#include <vector>
#include <numeric>

class Producer_Perforate : public Perforate{

 private:

    std::vector<double> _flow;

 public:

 Producer_Perforate() : Perforate(){};
    
 Producer_Perforate(const int& fluids_quantity) : Perforate(){
        _flow = std::vector<double>(fluids_quantity);
    };
    
    void flow(const int fluid_index, const double flow){_flow[fluid_index] = flow;};

    const std::vector<double>& flow() const {return _flow;};

    const double totalFlow() const override{
        double total_flow = 0;
        for(auto& fluid_flow : _flow) total_flow+=fluid_flow;
        return total_flow;//std::accumulate(_flow.begin(),_flow.end(),0);
    };


    const std::string type() const override {return typeid(Producer_Perforate).name();};
    const double flow(const int fluid_index) const {return _flow[fluid_index];};
};



#endif /* ṔRODUCER_PERFORATE_H */
