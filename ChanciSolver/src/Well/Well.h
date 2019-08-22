#ifndef WELL_H
#define WELL_H

#include <memory>
#include <vector>

#include "Value_Reader.h"
#include "Equation.h"
#include "Fluid.h"

#include "Mesh.h"
#include "Perforate.h"

#include "Operative_Condition.h"

class Well : public Equation<Well>{

 protected:

    using Perforates_t = std::vector<std::shared_ptr<Perforate>>;
    
    
    std::string _type;
    int _index;
    double _radius;
    int _number_of_perforates;
    double _borehole_depth;
    std::vector<std::shared_ptr<Perforate>> _perforates;

    std::vector<double> _borehole_pressure;
    std::vector<double> _flow;

    std::shared_ptr<Operative_Condition> _operative_condition;

    int _operative_status = 0; // 2 - Pending Change
                     // 1 - Changed
                     // 0 - Stable

 public:

    using Perforate_iterator = Perforates_t::iterator;
    using Perforate_const_iterator = Perforates_t::const_iterator;
    
 Well(const int index): _index(index){};

    virtual ~Well() = default;
    
    virtual void perforate(Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type){

        _perforates = std::vector<std::shared_ptr<Perforate>>();
        
        _type = type;

        _flow = std::vector<double>();

        _borehole_pressure=std::vector<double>();
        
        Value_Reader::myRead(std::string("Please insert the well radius "), _radius, std::string("Please insert a valid input"));
        Value_Reader::myRead(std::string("Please insert the number of perforations "), _number_of_perforates, std::string("Please insert a valid input"));

        Equation<Well>::_status = false;

        _operative_condition = std::make_shared<Operative_Condition>();

        _operative_status = 2;
        
    };

    virtual void perforateFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Fluid>>& characterized_fluids, const std::string& type){

        std::string element;
        
        _perforates = std::vector<std::shared_ptr<Perforate>>();
        
        _type = type;

        _flow = std::vector<double>();

        _borehole_pressure=std::vector<double>();
        
        while(well_reader >> element){
            std::transform(element.begin(), element.end(), element.begin(), ::toupper);
            if(element == "RADIUS"){
                well_reader >> _radius;
            }else if(element == "NUMBER_OF_PERFORATIONS"){
                well_reader >> _number_of_perforates;
                break;
            }
        };

        Equation<Well>::_status = false;
        
        _operative_condition = std::make_shared<Operative_Condition>();
        
        _operative_status = 2;
    };

    void flow(const int& term, const double flow){
        _flow[term] = flow;
    };
    
    void boreholePressure(const int& term, const double boreholePressure) {
        _borehole_pressure[term]=boreholePressure;
    };

    const int& index() const {return _index;};
    const double& radius() const {return _radius;};
    const double& boreholeDepth() const {return _borehole_depth;};
    const double& boreholePressure(const int& term) const { return _borehole_pressure[term];};
    
    const double& flow(const int& term) const { return _flow[term];};

    const int& numberOfPerforates() const {return _number_of_perforates;};

    template<typename PerforationType> inline void insertPerforationsFromFile(std::ifstream& well_reader, Mesh& mesh){
        std::string element;
        double skin;
        
        std::shared_ptr<Perforate> aux_perforate;

        std::vector<int> position = std::vector<int>(3);

        for(int perforate=0; perforate<_number_of_perforates; ++perforate){
            aux_perforate = std::make_shared<PerforationType>(PerforationType());
            while(well_reader >> element){
                std::transform(element.begin(), element.end(),element.begin(), ::toupper);                
                if(element == "POSITION"){
                    for(int axis=0; axis<3; ++axis){
                        well_reader >> position[axis];
                    };
                }else if(element == "SKIN_FACTOR"){
                    well_reader >> skin;
                    break;
                };
            };
            
            aux_perforate->position(position[0], position[1], position[2]);
            aux_perforate->index(mesh.listCell(position[0],position[1],position[2]));
            aux_perforate->localIndex(perforate);
            aux_perforate->skin(skin);
            
            if(perforate == 0){
                auto first_perforate_cell = mesh.cell(aux_perforate->index());
                _borehole_depth = first_perforate_cell->depth();
            };
            
        };
    };

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
            aux_perforate->skin(skin);
            
            if(perforate == 0){
                auto first_perforate_cell = mesh.cell(aux_perforate->index());
                _borehole_depth = first_perforate_cell->depth();
            };
        };
    };

    virtual void updateProperties(const int term){
        _borehole_pressure.push_back(_borehole_pressure[term-1]);
        _flow.push_back(_flow[term-1]);
    };

    std::shared_ptr<Operative_Condition>& operativeCondition(){
        return _operative_condition;
    };

    void establish(const int& term, std::string timestamp){

        std::string type;
        double value;
        double next_change;
        
        if((timestamp == "change" || timestamp == "") && _status == 2){
            
            Value_Reader::myRead(std::string("Please insert the type of operative condition (Pressure or Flow)"), type, std::string("Please insert a valid option"));

            std::transform(type.begin(), type.end(),type.begin(), ::toupper);
            
            while(type != "PRESSURE" && type != "FLOW" && type != "SHUT"){
                Value_Reader::myRead(std::string("Please insert the type of operative condition (Pressure or Flow)"), type, std::string("Please insert a valid option"));
            };
            
            Value_Reader::myRead(std::string("Please insert the value of operative condition "), value, std::string("Please insert a valid option"));

            Value_Reader::myRead(std::string("Please insert the next time of change of operative condition "), next_change, std::string("Please insert a valid option"));

            _operative_condition->type(type);
            _operative_condition->value(value);
            _operative_condition->nextChange(next_change);

            if(type == "PRESSURE" || type == "SHUT"){
                Equation<Well>::_status = false;
            }else{
                Equation<Well>::_status = true;
            };

            _operative_status = 1; //
        };
            
    };

    void establishFromFile(std::ifstream& condition_reader, const int& term, std::string timestamp){

        std::string element;
        std::string type;
        double value;
        double next_change;
        
        if((timestamp == "change" || timestamp == "") && _status == 2){

            while(condition_reader >> element){

                std::transform(element.begin(), element.end(),element.begin(), ::toupper);

                if(element == "TYPE"){
                    condition_reader >> type;
                    std::transform(type.begin(), type.end(),type.begin(), ::toupper);
                }else if(element == "VALUE"){
                    condition_reader >> value;
                }else if(element == "NEXT_TIME"){
                    condition_reader >> next_change;
                    break;
                };
            };

            _operative_condition->type(type);
            _operative_condition->value(value);
            _operative_condition->nextChange(next_change);
            
            if(type == "PRESSURE"  || type == "SHUT"){
                Equation<Well>::_status = false;
            }else{
                Equation<Well>::_status = true;
            };

            _operative_status = 1;
        };
            
    };

    void operativeStatus(const int operative_status){_operative_status=operative_status;};
    const int& operativeStatus() const {return _operative_status;};
    
    Perforate_iterator begin() {return _perforates.begin();};
    Perforate_iterator end()   {return _perforates.end();};

    Perforate_const_iterator begin()  const {return _perforates.begin();};
    Perforate_const_iterator end()    const {return _perforates.end();};
    Perforate_const_iterator cbegin() const {return _perforates.cbegin();};
    Perforate_const_iterator cend()   const {return _perforates.cend();};
    
};



#endif /* WELL_H */
