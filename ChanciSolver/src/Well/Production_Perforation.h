#ifndef PRODUCTION_PERFORATION_H
#define PRODUCTION_PERFORATION_H

#include "Perforation.h"
#include <vector>
#include <numeric>

class Production_Perforation : public Perforation{

 private:

    std::vector<double> _flow;

 public:

 Production_Perforation() : Perforation(){};
    
 Production_Perforation(const int& phases_quantity) : Perforation(){
        _flow = std::vector<double>(phases_quantity);
    };
    
    void flow(const int phase_index, const double flow){_flow[phase_index] = flow;};

    const std::vector<double>& flow() const {return _flow;};

    const double totalFlow() const override{
        double total_flow = 0;
        for(auto& phase_flow : _flow) total_flow+=phase_flow;
        return total_flow;//std::accumulate(_flow.begin(),_flow.end(),0);
    };


    const std::string type() const override {return typeid(Production_Perforation).name();};
    const double flow(const int phase_index) const {return _flow[phase_index];};
};



#endif /* PRODUCTION_PERFORATION_H */
