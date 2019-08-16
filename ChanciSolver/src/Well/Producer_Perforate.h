#ifndef PRODUCER_PERFORATE_H
#define PRODUCER_PERFORATE_H

#include "Perforate.h"
#include <vector>
#include <algorithm>

class Producer_Perforate : public Perforate{

 private:

    std::vector<int> _ignore_indexes;
    std::vector<double> _flow;

 public:

 Producer_Perforate() : Perforate(){};
    
    void flow(const int fluid_index, const double flow){_flow[fluid_index] = flow;};

    const std::vector<double>& flow() const {return _flow;};

    const double totalFlow() const override{
        double total=0;
        for(int fluid_index=0; fluid_index<_flow.size(); ++fluid_index){
            auto result = std::find(std::begin(_ignore_indexes), std::end(_ignore_indexes), fluid_index);
            if(result != std::end(_ignore_indexes)){
                total+=_flow[fluid_index];
            };
        };
    };


    const std::string type() const override {return typeid(Producer_Perforate).name();};
};



#endif /* ṔRODUCER_PERFORATE_H */
