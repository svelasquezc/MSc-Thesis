#ifndef WELL_H
#define WELL_H

#include <memory>
#include <vector>

#include "Value_Reader.h"
#include "Equation.h"

class Well : std::enable_shared_from_this<Well>{

 private:
    
    std::string _type;
    int _index;
    double _radius;
    int _number_of_perforates;
    double _borehole_depth;
    std::vector<Perforate> _perforates;

    std::vector<double> _borehole_pressure;
    std::vector<double> _flow;

    Equation<Well> _equation;

 public:
    
    void perforate(std::vector<std::shared_ptr<Fluid>>& characterized_fluids){
        Value_Reader::myRead(std::string("Please insert the type of well "), _type, std::string("Please insert a valid input"));

        Value_Reader::myRead(std::string("Please insert the well radius "), _radius, std::string("Please insert a valid input"));
        Value_Reader::myRead(std::string("Please insert the number of perforates "), _number_of_perforates, std::string("Please insert a valid input"));

        for(int perforate_index=0; perforate_index < _number_of_perforates; ++perforate_index){
            
        };
    };
};



#endif /* WELL_H */
