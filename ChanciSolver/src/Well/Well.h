#ifndef WELL_H
#define WELL_H

#include <memory>
#include <vector>

#include "Value_Reader.h"
#include "Equation.h"
#include "Fluid.h"

#include "Perforate.h"
#include "Producer_Perforate.h"
#include "Injector_Perforate.h"

class Well : public Equation<Well>, public std::enable_shared_from_this<Well>{

 private:
    
    std::string _type;
    int _index;
    double _radius;
    int _number_of_perforates;
    double _borehole_depth;
    std::vector<std::unique_ptr<Perforate>> _perforates;

    std::vector<double> _borehole_pressure;
    std::vector<double> _flow;

    

 public:
 Well() : Equation<Well>(shared_from_this()){};
    virtual void perforate(std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type){
        _type = type;
        
        Value_Reader::myRead(std::string("Please insert the well radius "), _radius, std::string("Please insert a valid input"));
        Value_Reader::myRead(std::string("Please insert the number of perforates "), _number_of_perforates, std::string("Please insert a valid input"));

        for(int perforate_index=0; perforate_index < _number_of_perforates; ++perforate_index){
            
        };
    };

    void boreholePressure(const int& term, const double boreholePressure) {
        _borehole_pressure[term]=boreholePressure;
    };
    
    const double& boreholePressure(const int& term) const { return _borehole_pressure[term];};

    void flow(const int& term, const double flow) {
        _flow[term]=flow;
    };
    
    const double& flow(const int& term) const { return _flow[term];};

    virtual ~Well() = default;
    
};



#endif /* WELL_H */
