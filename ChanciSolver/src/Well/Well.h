#ifndef WELL_H
#define WELL_H

#include <memory>
#include <vector>

#include "Value_Reader.h"
#include "Equation.h"
#include "Fluid.h"

#include "Mesh.h"
#include "Perforate.h"

class Well : public Equation<Well>{

 protected:
    
    std::string _type;
    int _index;
    double _radius;
    int _number_of_perforates;
    double _borehole_depth;
    std::vector<std::shared_ptr<Perforate>> _perforates;

    std::vector<double> _borehole_pressure;
    std::vector<double> _flow;

    

 public:
    Well(){};

    virtual ~Well() = default;
    
    virtual void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type){

        _perforates = std::vector<std::shared_ptr<Perforate>>();
        
        _type = type;

        _flow = std::vector<double>();

        _borehole_pressure=std::vector<double>();
        
        Value_Reader::myRead(std::string("Please insert the well radius "), _radius, std::string("Please insert a valid input"));
        Value_Reader::myRead(std::string("Please insert the number of perforates "), _number_of_perforates, std::string("Please insert a valid input"));

    };

    void boreholePressure(const int& term, const double boreholePressure) {
        _borehole_pressure[term]=boreholePressure;
    };

    const int& index() const {return _index;};
    const double& radius() const {return _radius;};
    const double& boreholeDepth() const {return _borehole_depth;};
    const double& boreholePressure(const int& term) const { return _borehole_pressure[term];};

    void flow(const int& term, const double flow) {
        _flow[term]=flow;
    };
    
    const double& flow(const int& term) const { return _flow[term];};

    const int& numberOfPerforates() const {return _number_of_perforates;};

    template<typename PerforationType> inline void insertPerforations(Mesh& mesh){

        double skin;
        
        std::ostringstream ss = std::ostringstream();
        const std::string axisnames[3]={"x", "y", "z"};
        
        std::shared_ptr<Perforate> aux_perforate;

        std::vector<int> position = std::vector<int>(3);

        for(int perforate=0; perforate<_number_of_perforates; ++perforate){
            aux_perforate = std::make_shared<PerforationType>(PerforationType());
            for(int axis=0; axis<3; ++axis){
                
                ss << "Please insert perforate "<< perforate + 1 << " position in axis " << axisnames[axis] << ": ";
                Value_Reader::myRead(ss.str(), position[axis], std::string("Please insert a valid option"));
                ss.str("");
                ss.clear();
            };
            
            ss << "Please insert perforate "<< perforate + 1 << " skin factor: ";        
            Value_Reader::myRead(ss.str(), skin, std::string("Please insert a valid option"));
            
            aux_perforate->position(position[0], position[1], position[2]);
            aux_perforate->index(mesh.listCell(position[0],position[1],position[2]));
            aux_perforate->localIndex(perforate);
            
            
        };
    };
    
};



#endif /* WELL_H */
