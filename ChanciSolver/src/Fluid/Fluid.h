#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>


#include "Measured_Property.h"
#include "Equation.h"

class Fluid : public Equation<Fluid>{
    
 private:
    
    int _index;
    static int _count_of_principals;
    static int _count_of_fluids    ;
    std::string _type;
    std::map<const std::string,std::vector<double>> _pressure;
    std::map<const std::string,std::vector<double>> _density;
    std::map<const std::string,std::vector<double>> _saturation;
    std::map<const std::string,std::vector<double>> _viscosity;
    std::map<const std::string,std::vector<double>> _volumetric_factor;
    std::map<const std::string,std::vector<double>> _potential;
    std::map<const std::string,std::vector<double>> _relative_permeability;
    double _standard_conditions_density;
    std::unique_ptr<Measured_Property> _measured_volumetric_factor;
    std::unique_ptr<Measured_Property> _measured_viscosity;
    mutable bool _principal=false;

 public:
    Fluid(){};
 Fluid(const int index) : _index(index){};
    void characterize(int& cells_number);
    void characterizeFromFile(std::ifstream& fluid_reader, int& cells_number);
    void updateProperties(const std::string& term);
    //void calculate(int& term, int& _cellindex);
    const std::string& print() const;
    const int& index() const;
    //sets
    void volumetricFactor(const std::string& _term, const int& _cellindex, const double pressure);
    void viscosity(const std::string& _term, const int& _cellindex, const double pressure);
    void potential(const std::string& _term, const int& _cellindex, const double _gravity, const double& depth);

    void density(const std::string& term, const int& cell_index, const double density){_density[term][cell_index]=density;}
    void pressure(const std::string& term, const int& cell_index, const double pressure){_pressure[term][cell_index]=pressure;}

    void saturation(const std::string& term, const int& cell_index, const double saturation){_saturation[term][cell_index]=saturation;}

    void relativePermeability(const std::string& term, const int& cell_index, const double relative_permeability){_relative_permeability[term][cell_index]=relative_permeability;}
    
    //gets
    const double& standardConditionsDensity() const {return _standard_conditions_density;};
    const double& pressure             (const std::string& _term, const int _cell_index) const {return _pressure.find(_term)->second[_cell_index];};             
    const double& density              (const std::string& _term, const int _cell_index) const {return _density.find(_term)->second[_cell_index];};              
    const double& saturation           (const std::string& _term, const int _cell_index) const {return _saturation.find(_term)->second[_cell_index];};           
    const double& viscosity            (const std::string& _term, const int _cell_index) const {return _viscosity.find(_term)->second[_cell_index];};            
    const double& volumetricFactor     (const std::string& _term, const int _cell_index) const {return _volumetric_factor.find(_term)->second[_cell_index];};    
    const double& potential            (const std::string& _term, const int _cell_index) const {return _potential.find(_term)->second[_cell_index];};            
    const double& relativePermeability (const std::string& _term, const int _cell_index) const {return _relative_permeability.find(_term)->second[_cell_index];};
    const bool  & principal            ()                                       const {return _principal;};
    
    static const int& countOfPrincipals(){return _count_of_principals;};
    static const int& countOfFluids    (){return _count_of_fluids;};
    
};

void Fluid::characterize(int& cells_number){

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();
    std::ostringstream ss = std::ostringstream();
    
    _pressure              ["N"] = std::vector<double>(cells_number);
    _density               ["N"] = std::vector<double>(cells_number);
    _saturation            ["N"] = std::vector<double>(cells_number);
    _viscosity             ["N"] = std::vector<double>(cells_number);
    _volumetric_factor     ["N"] = std::vector<double>(cells_number);
    _relative_permeability ["N"] = std::vector<double>(cells_number);
    _potential             ["N"] = std::vector<double>(cells_number);

    _pressure              ["K"] = std::vector<double>(cells_number);
    _density               ["K"] = std::vector<double>(cells_number);
    _saturation            ["K"] = std::vector<double>(cells_number);
    _viscosity             ["K"] = std::vector<double>(cells_number);
    _volumetric_factor     ["K"] = std::vector<double>(cells_number);
    _relative_permeability ["K"] = std::vector<double>(cells_number);
    _potential             ["K"] = std::vector<double>(cells_number);

    Value_Reader::myRead(std::string("Please insert the type of fluid "), _type, std::string("Please insert a valid input"));
    
    // Reading of measured properties (PVT Table)
    ref_name  << _type << " Pressure";
    ref_value << _type << " Volumetric Factor";
    _measured_volumetric_factor = std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
    ref_value.str("");
    ref_value.clear();
    
    ref_value << _type << " Viscosity";
    _measured_viscosity = std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
    
    _measured_volumetric_factor->readMe();
    _measured_viscosity->readMe();
    
    // Reading of standard conditions density
    Value_Reader::myRead(std::string("Please insert the standard conditions density for the fluid "), _standard_conditions_density, std::string("Please insert a valid input"));
    
    // Initial Conditions for the fluid
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        ss << "Please insert initial pressure for the "<< cellindex+1 << " cell [Pa]";
        Value_Reader::myRead(ss.str(), _pressure["N"][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
    };
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        ss << "Please insert initial saturation for the "<< cellindex+1 << " cell [-]";
        Value_Reader::myRead(ss.str(), _saturation["N"][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
    };

    if(_count_of_principals < 1) {
        bool aux_principal=false;
        Value_Reader::myRead(std::string("Please insert 1 or 0 for principal or not principal fluid (bool) "), aux_principal, std::string("Please insert a valid input"));
        _principal=aux_principal;
        if(_principal){
            ++_count_of_principals;
        }
    }else{
        _principal=false;
    };

    _index = _count_of_fluids;
    ++_count_of_fluids;

    Equation<Fluid>::_status = true;
};

void Fluid::characterizeFromFile(std::ifstream& fluid_reader, int& cells_number){
    std::string element;

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();

    _pressure              ["N"] = std::vector<double>(cells_number);
    _density               ["N"] = std::vector<double>(cells_number);
    _saturation            ["N"] = std::vector<double>(cells_number);
    _viscosity             ["N"] = std::vector<double>(cells_number);
    _volumetric_factor     ["N"] = std::vector<double>(cells_number);
    _relative_permeability ["N"] = std::vector<double>(cells_number);
    _potential             ["N"] = std::vector<double>(cells_number);

    _pressure              ["K"] = std::vector<double>(cells_number);
    _density               ["K"] = std::vector<double>(cells_number);
    _saturation            ["K"] = std::vector<double>(cells_number);
    _viscosity             ["K"] = std::vector<double>(cells_number);
    _volumetric_factor     ["K"] = std::vector<double>(cells_number);
    _relative_permeability ["K"] = std::vector<double>(cells_number);
    _potential             ["K"] = std::vector<double>(cells_number);
    
    // Reading of measured properties (PVT Table)
    ref_name  << _type << " Pressure";
    ref_value << _type << " Volumetric Factor";
    _measured_volumetric_factor = std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
    ref_value.str("");
    ref_value.clear();
    
    ref_value << _type << " Viscosity";
    _measured_viscosity = std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
    
    
    while(fluid_reader>>element){
        
        std::transform(element.begin(), element.end(),element.begin(), ::toupper);
        
        if(element == "TYPE"){
            
            fluid_reader>>_type;
            
        }else if(element == "VOLUME_FACTOR"){
            
            _measured_volumetric_factor->readFromFile(fluid_reader);
            
        }else if(element == "VISCOSITY"){
            
            _measured_viscosity->readFromFile(fluid_reader);
            
        }else if(element == "STANDARD_CONDITIONS_DENSITY"){
            
            fluid_reader>>_standard_conditions_density;
            
        }else if(element == "INITIAL_PRESSURE"){
            
            for(int cellindex=0; cellindex<cells_number; ++cellindex){
                fluid_reader>>_pressure["N"][cellindex];
            };
            
        }else if(element == "INITIAL_SATURATION"){
            
            for(int cellindex=0; cellindex<cells_number; ++cellindex){
                fluid_reader>>_saturation["N"][cellindex];
            };
            
        }else if(element == "PRINCIPAL"){
            
            fluid_reader>>_principal;
            if(_count_of_principals < 1) {
                if(_principal){
                    ++_count_of_principals;
                }
            }else{
                _principal=false;
            };
            
            break;
        }
    };
    Equation<Fluid>::_status = true;
};

void Fluid::updateProperties(const std::string& term){
    if(term == "K"){
        _pressure["K"]             =_pressure["N"];
        _potential["K"]            =_potential["N"];
        _density["K"]              =_density["N"];
        _saturation["K"]           =_saturation["N"];
        _viscosity["K"]            =_viscosity["N"];
        _volumetric_factor["K"]    =_volumetric_factor["N"];
        _relative_permeability["K"]=_relative_permeability["N"];
    }else{
        _pressure["N"]             =_pressure["K"];
        _potential["N"]            =_potential["K"];
        _density["N"]              =_density["K"];
        _saturation["N"]           =_saturation["K"];
        _viscosity["N"]            =_viscosity["K"];
        _volumetric_factor["N"]    =_volumetric_factor["K"];
        _relative_permeability["N"]=_relative_permeability["K"];
    }
};

void Fluid::volumetricFactor(const std::string& _term, const int& _cellindex, const double pressure){
    /*
      Here we need to calculate the restrictions to the flow equations (Volume restriction and Capilarity)

      //Saturation calculation(Volume Restriction)

      //Pressure calculation (Capilarity)

      */
    //Properties Calculation
    _volumetric_factor[_term][_cellindex] =
        _measured_volumetric_factor->interpolate(_pressure[_term][_cellindex]);
};

void Fluid::viscosity(const std::string& _term, const int& _cellindex, const double pressure){
    _viscosity[_term][_cellindex] =
        _measured_viscosity->interpolate(pressure);
};

const std::string& Fluid::print() const{
    return _type;
};

const int& Fluid::index() const{
    return _index;
};

void Fluid::potential(const std::string& _term, const int& _cellindex, const double gravity, const double& depth){
    _potential[_term][_cellindex] =
        _pressure[_term][_cellindex] + _density[_term][_cellindex]*gravity*depth;
};

//void {density[_term][_cellindex] = standard_conditions_density;};

#endif /* FLUID_H */
