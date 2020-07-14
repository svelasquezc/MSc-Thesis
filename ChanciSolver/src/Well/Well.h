#ifndef WELL_H
#define WELL_H

#include <memory>
#include <vector>

#include "Value_Reader.h"
#include "Equation.h"
#include "Phase.h"

#include "Mesh.h"
#include "Perforation.h"

#include "Operative_Condition.h"

class Well : public Equation<Well>{

protected:

    using Perforations_t = std::vector<std::shared_ptr<Perforation>>;
    
    std::string _type;
    int _index;
    double _radius;
    int _perforations_quantity;
    double _borehole_depth;
    double _average_density;
    std::vector<std::shared_ptr<Perforation>> _perforations;

    std::vector<double> _borehole_pressure;
    std::vector<double> _flow;

    std::shared_ptr<Operative_Condition> _operative_condition;

    int _operative_status = 0; // 2 - Pending Change
    // 1 - Changed
    // 0 - Stable

public:

    using Perforation_iterator = Perforations_t::iterator;
    using Perforation_const_iterator = Perforations_t::const_iterator;
    
    Well(const int index): _index(index){};

    virtual ~Well() = default;
    
    virtual void perforation(Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type){

        _perforations = std::vector<std::shared_ptr<Perforation>>();
        
        _type = type;

        _flow = std::vector<double>(2,0.0);
        
        _borehole_pressure=std::vector<double>(2,0.0);
        
        Value_Reader::myRead(std::string("Please insert the well radius "), _radius, std::string("Please insert a valid input"));
        Value_Reader::myRead(std::string("Please insert the number of perforations "), _perforations_quantity, std::string("Please insert a valid input"));

        Equation<Well>::_status = false;

        _operative_condition = std::make_shared<Operative_Condition>();

        _operative_status = 2;
        
    };

    virtual void perforationFromFile(std::ifstream& well_reader, Mesh& mesh, std::vector<std::shared_ptr<Phase>>& characterized_phases, const std::string& type){

        std::string element;
        
        _perforations = std::vector<std::shared_ptr<Perforation>>();
        
        _type = type;
        
        _flow = std::vector<double>(2,0.0);
        
        _borehole_pressure=std::vector<double>(2,0.0);
        
        while(well_reader >> element){
            std::transform(element.begin(), element.end(), element.begin(), ::toupper);
            if(element == "RADIUS"){
                well_reader >> _radius;
            }else if(element == "PERFORATIONS_QUANTITY"){
                well_reader >> _perforations_quantity;
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

    const int& numberOfPerforations() const {return _perforations_quantity;};

    template<typename PerforationType> inline void insertPerforationsFromFile(std::ifstream& well_reader, Mesh& mesh, const int& phases_quantity){
        std::string element;
        double skin;
        
        std::shared_ptr<Perforation> aux_perforation;

        std::vector<int> location = std::vector<int>(3);

        for(int perforation=0; perforation<_perforations_quantity; ++perforation){
            aux_perforation = std::make_shared<PerforationType>(PerforationType(phases_quantity));
            while(well_reader >> element){
                std::transform(element.begin(), element.end(),element.begin(), ::toupper);                
                if(element == "LOCATION"){
                    for(int axis=0; axis<3; ++axis){
                        well_reader >> location[axis];
                        --location[axis];
                    };
                }else if(element == "SKIN_FACTOR"){
                    well_reader >> skin;
                    break;
                };
            };
            
            aux_perforation->location(location[0], location[1], location[2]);
            aux_perforation->index(mesh.listCell(location[0],location[1],location[2]));
            aux_perforation->localIndex(perforation);
            aux_perforation->skin(skin);
            
            if(perforation == 0){
                auto first_perforation_cell = mesh.cell(aux_perforation->index());
                _borehole_depth = first_perforation_cell->depth();
            };

            _perforations.push_back(aux_perforation);
            
        };
    };

    template<typename PerforationType> inline void insertPerforations(Mesh& mesh, const int& phases_quantity){

        double skin;
        
        std::ostringstream ss = std::ostringstream();
        const std::string axisnames[3]={"x", "y", "z"};
        
        std::shared_ptr<Perforation> aux_perforation;

        std::vector<int> location = std::vector<int>(3);

        for(int perforation=0; perforation<_perforations_quantity; ++perforation){
            aux_perforation = std::make_shared<PerforationType>(PerforationType(phases_quantity));
            for(int axis=0; axis<3; ++axis){
                
                ss << "Please insert perforation "<< perforation + 1 << " location in axis " << axisnames[axis] << ": ";
                Value_Reader::myRead(ss.str(), location[axis], std::string("Please insert a valid option"));
                ss.str("");
                ss.clear();
            };
            
            ss << "Please insert perforation "<< perforation + 1 << " skin factor: ";        
            Value_Reader::myRead(ss.str(), skin, std::string("Please insert a valid option"));
            
            aux_perforation->location(location[0], location[1], location[2]);
            aux_perforation->index(mesh.listCell(location[0],location[1],location[2]));
            aux_perforation->localIndex(perforation);
            aux_perforation->skin(skin);
            
            if(perforation == 0){
                auto first_perforation_cell = mesh.cell(aux_perforation->index());
                _borehole_depth = first_perforation_cell->depth();
            };

            _perforations.push_back(aux_perforation);
        };
    };

    virtual void updateProperties(const int& term){
        if(term==1){
            _borehole_pressure[1]=_borehole_pressure[0];
            _flow[1]=_flow[0];
        }else{
            _borehole_pressure[0]=_borehole_pressure[1];
            _flow[0]=_flow[1];
        }
        
        for(auto perforation : _perforations){
            perforation->updateProperties(term);
        }
    };

    std::shared_ptr<Operative_Condition>& operativeCondition(){
        return _operative_condition;
    };

    void establish(const int& term, std::string timestamp){

        std::string type;
        double value;
        double next_change;
        
        if((timestamp == "change" || timestamp == "") && _operative_status == 2){
            
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
                if(type == "PRESSURE"){
                    _borehole_pressure[term] = value;
                };
            }else{
                Equation<Well>::_status = true;
                _flow[term] = value;
            };

            _operative_status = 1; //
        };
            
    };

    void establishFromFile(std::ifstream& condition_reader, const int& term, std::string timestamp){

        std::string element;
        std::string type;
        double value;
        double next_change;
        
        if((timestamp == "change" || timestamp == "") && _operative_status == 2){

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

            if(type == "PRESSURE" && type == "SHUT" && type != "FLOW"){
                throw std::invalid_argument("Type of operative condition must be FLOW, or PRESSURE or SHUT");
            };
            
            _operative_condition->type(type);
            _operative_condition->value(value);
            _operative_condition->nextChange(next_change);
            
            if(type == "PRESSURE" || type == "SHUT"){
                Equation<Well>::_status = false;
                if(type == "PRESSURE"){
                    _borehole_pressure[term] = value;
                };
            }else{
                Equation<Well>::_status = true;
                _flow[term] = value;
            };

            _operative_status = 1;
        };
            
    };

    virtual const std::string type() const {return typeid(Well).name();};

    void operativeStatus(const int operative_status){_operative_status=operative_status;};
    const int& operativeStatus() const {return _operative_status;};
    
    Perforation_iterator begin() {return _perforations.begin();};
    Perforation_iterator end()   {return _perforations.end();};

    Perforation_const_iterator begin()  const {return _perforations.begin();};
    Perforation_const_iterator end()    const {return _perforations.end();};
    Perforation_const_iterator cbegin() const {return _perforations.cbegin();};
    Perforation_const_iterator cend()   const {return _perforations.cend();};
};



#endif /* WELL_H */
