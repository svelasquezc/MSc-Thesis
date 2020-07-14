#ifndef PHASE_H
#define PHASE_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>


#include "Property_Measure.h"
#include "Equation.h"

class Phase : public Equation<Phase>{
    
private:
    
    int _index;
    static int _main_phase_counter;
    static int _count_of_phases    ;
    std::string _type;
    std::vector<std::vector<double>> _pressure;
    std::vector<std::vector<double>> _density;
    std::vector<std::vector<double>> _saturation;
    std::vector<std::vector<double>> _viscosity;
    std::vector<std::vector<double>> _formation_volume_factor;
    std::vector<std::vector<double>> _potential;
    std::vector<std::vector<double>> _relative_permeability;
    double _standard_conditions_density;
    std::unique_ptr<Property_Measure> _formation_volume_factor_measure;
    std::unique_ptr<Property_Measure> _viscosity_measure;
    mutable bool _main=false;

public:
    Phase(){};
    Phase(const int index) : _index(index){};
    void characterize(int& cells_quantity);
    void characterizeFromFile(std::ifstream& phase_reader, int& cells_quantity);
    void updateProperties(const int& term);
    //void calculate(int& term, int& _cellindex);
    const std::string& print() const;
    const int& index() const;
    //sets
    void formationVolumeFactor(const int& _term, const int& _cellindex, const double pressure);
    void viscosity(const int& _term, const int& _cellindex, const double pressure);
    void potential(const int& _term, const int& _cellindex, const double _gravity, const double& depth);

    void density(const int& term, const int& cell_index, const double density){_density[term][cell_index]=density;}
    void pressure(const int& term, const int& cell_index, const double pressure){_pressure[term][cell_index]=pressure;}

    void saturation(const int& term, const int& cell_index, const double saturation){_saturation[term][cell_index]=saturation;}

    void relativePermeability(const int& term, const int& cell_index, const double relative_permeability){_relative_permeability[term][cell_index]=relative_permeability;}
    
    //gets
    const double& standardConditionsDensity() const {return _standard_conditions_density;};
    const double& pressure             (const int& _term, const int _cell_index) const {return _pressure[_term][_cell_index];};             
    const double& density              (const int& _term, const int _cell_index) const {return _density[_term][_cell_index];};              
    const double& saturation           (const int& _term, const int _cell_index) const {return _saturation[_term][_cell_index];};           
    const double& viscosity            (const int& _term, const int _cell_index) const {return _viscosity[_term][_cell_index];};            
    const double& formationVolumeFactor     (const int& _term, const int _cell_index) const {return _formation_volume_factor[_term][_cell_index];};    
    const double& potential            (const int& _term, const int _cell_index) const {return _potential[_term][_cell_index];};            
    const double& relativePermeability (const int& _term, const int _cell_index) const {return _relative_permeability[_term][_cell_index];};
    const bool  & main            ()                                       const {return _main;};

    const std::vector<double>& pressure (const int& term) const {return _pressure[term];};
    const std::vector<double>& density (const int& term) const {return _density[term];};
    const std::vector<double>& saturation (const int& term) const {return _saturation[term];};
    const std::vector<double>& viscosity (const int& term) const {return _viscosity[term];};
    const std::vector<double>& formationVolumeFactor (const int& term) const {return _formation_volume_factor[term];};
    const std::vector<double>& potential (const int& term) const {return _potential[term];};
    const std::vector<double>& relativePermeability (const int& term) const {return _relative_permeability[term];};
    
    static const int& countOfMains(){return _main_phase_counter;};
    static const int& countOfPhases    (){return _count_of_phases;};


    template<typename VTKType> void insert(VTKType& vtkcontext, const int term){

        std::ostringstream ss_pressure;
        std::ostringstream ss_saturation;
        std::ostringstream ss_viscosity;
        std::ostringstream ss_density;
        std::ostringstream ss_potential;
        std::ostringstream ss_formation_volume_factor;
        std::ostringstream ss_relative_permeability;

        ss_pressure << _type << "_Pressure";
        ss_saturation << _type << "_Saturation";
        ss_viscosity << _type << "_Viscosity";
        ss_density << _type << "_Density";
        ss_potential << _type << "_Potential";
        ss_formation_volume_factor << _type << "_Formation_Volume_Factor";
        ss_relative_permeability << _type << "_Relative_Permeability";
        
        vtkcontext.appendScalar(ss_pressure.str(), _pressure[term]);
        vtkcontext.appendScalar(ss_saturation.str(), _saturation[term]);
        vtkcontext.appendScalar(ss_density.str(), _density[term]);
        vtkcontext.appendScalar(ss_formation_volume_factor.str(), _formation_volume_factor[term]);
        vtkcontext.appendScalar(ss_viscosity.str(), _viscosity[term]);
        vtkcontext.appendScalar(ss_potential.str(), _potential[term]);
        vtkcontext.appendScalar(ss_relative_permeability.str(), _relative_permeability[term]);
    };
    
};

void Phase::characterize(int& cells_quantity){

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();
    std::ostringstream ss = std::ostringstream();
    
    _pressure              = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _density               = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _saturation            = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _viscosity             = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _formation_volume_factor     = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _relative_permeability = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _potential             = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));

    Value_Reader::myRead(std::string("Please insert the type of phase "), _type, std::string("Please insert a valid input"));
    
    // Reading of measured properties (PVT Table)
    ref_name  << _type << " Pressure";
    ref_value << _type << " Volumetric Factor";
    _formation_volume_factor_measure = std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
    ref_value.str("");
    ref_value.clear();
    
    ref_value << _type << " Viscosity";
    _viscosity_measure = std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
    
    _formation_volume_factor_measure->readMe();
    _viscosity_measure->readMe();
    
    // Reading of standard conditions density
    Value_Reader::myRead(std::string("Please insert the standard conditions density for the phase "), _standard_conditions_density, std::string("Please insert a valid input"));
    
    // Initial Conditions for the phase
    for(int cellindex=0; cellindex<cells_quantity; ++cellindex){
        ss << "Please insert initial pressure for the "<< cellindex+1 << " cell [Pa]";
        Value_Reader::myRead(ss.str(), _pressure[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
    };
    for(int cellindex=0; cellindex<cells_quantity; ++cellindex){
        ss << "Please insert initial saturation for the "<< cellindex+1 << " cell [-]";
        Value_Reader::myRead(ss.str(), _saturation[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
    };

    if(_main_phase_counter < 1) {
        bool aux_main=false;
        Value_Reader::myRead(std::string("Please insert 1 or 0 for main or not main phase (bool) "), aux_main, std::string("Please insert a valid input"));
        _main=aux_main;
        if(_main){
            ++_main_phase_counter;
        }
    }else{
        _main=false;
    };

    _index = _count_of_phases;
    ++_count_of_phases;

    Equation<Phase>::_status = true;
};

void Phase::characterizeFromFile(std::ifstream& phase_reader, int& cells_quantity){
    std::string element;

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();

    _pressure              = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _density               = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _saturation            = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _viscosity             = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _formation_volume_factor     = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _relative_permeability = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    _potential             = std::vector<std::vector<double>>(2,std::vector<double>(cells_quantity));
    
    // Reading of measured properties (PVT Table)
    ref_name  << _type << " Pressure";
    ref_value << _type << " Volumetric Factor";
    _formation_volume_factor_measure = std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
    ref_value.str("");
    ref_value.clear();
    
    ref_value << _type << " Viscosity";
    _viscosity_measure = std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
    
    
    while(phase_reader>>element){
        
        std::transform(element.begin(), element.end(),element.begin(), ::toupper);
        
        if(element == "TYPE"){
            
            phase_reader>>_type;
            
        }else if(element == "FORMATION_VOLUME_FACTOR"){
            
            _formation_volume_factor_measure->readFromFile(phase_reader);
            
        }else if(element == "VISCOSITY"){
            
            _viscosity_measure->readFromFile(phase_reader);
            
        }else if(element == "STANDARD_CONDITIONS_DENSITY"){
            
            phase_reader>>_standard_conditions_density;
            
        }else if(element == "INITIAL_PRESSURE"){
            
            for(int cellindex=0; cellindex<cells_quantity; ++cellindex){
                phase_reader>>_pressure[0][cellindex];
            };
            
        }else if(element == "INITIAL_SATURATION"){
            
            for(int cellindex=0; cellindex<cells_quantity; ++cellindex){
                phase_reader>>_saturation[0][cellindex];
            };
            
        }else if(element == "MAIN"){
            
            phase_reader>>_main;
            if(_main_phase_counter < 1) {
                if(_main){
                    ++_main_phase_counter;
                }
            }else{
                _main=false;
            };
            
            break;
        }
    };
    Equation<Phase>::_status = true;
};

void Phase::updateProperties(const int& term){
    if(term == 1){
        _pressure[1]             =_pressure[0];
        _potential[1]            =_potential[0];
        _density[1]              =_density[0];
        _saturation[1]           =_saturation[0];
        _viscosity[1]            =_viscosity[0];
        _formation_volume_factor[1]    =_formation_volume_factor[0];
        _relative_permeability[1]=_relative_permeability[0];
    }else{
        _pressure[0]             =_pressure[1];
        _potential[0]            =_potential[1];
        _density[0]              =_density[1];
        _saturation[0]           =_saturation[1];
        _viscosity[0]            =_viscosity[1];
        _formation_volume_factor[0]    =_formation_volume_factor[1];
        _relative_permeability[0]=_relative_permeability[1];
    }
};

void Phase::formationVolumeFactor(const int& _term, const int& _cellindex, const double pressure){
    /*
      Here we need to calculate the restrictions to the flow equations (Volume restriction and Capilarity)

      //Saturation calculation(Volume Restriction)

      //Pressure calculation (Capilarity)

      */
    //Properties Calculation
    _formation_volume_factor[_term][_cellindex] =
        _formation_volume_factor_measure->interpolate(_pressure[_term][_cellindex]);
};

void Phase::viscosity(const int& _term, const int& _cellindex, const double pressure){
    _viscosity[_term][_cellindex] =
        _viscosity_measure->interpolate(pressure);
};

const std::string& Phase::print() const{
    return _type;
};

const int& Phase::index() const{
    return _index;
};

void Phase::potential(const int& _term, const int& _cellindex, const double gravity, const double& depth){
    _potential[_term][_cellindex] =
        _pressure[_term][_cellindex] + _density[_term][_cellindex]*gravity*depth;
};

//void {density[_term][_cellindex] = standard_conditions_density;};

#endif /* PHASE_H */
