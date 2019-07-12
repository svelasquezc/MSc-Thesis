#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <memory>
#include <string>
#include <sstream>

#include "Measured_Property.h"

class Fluid : protected Value_Reader{
    
 private:
    
    int _index;
    std::string _type;
    std::vector<std::vector<double>> _pressure;
    std::vector<std::vector<double>> _density;
    std::vector<std::vector<double>> _saturation;
    std::vector<std::vector<double>> _viscosity;
    std::vector<std::vector<double>> _volumetric_factor;
    std::vector<std::vector<double>> _potential;
    std::vector<std::vector<double>> _relative_permeability;
    double _standard_conditions_density;
    std::unique_ptr<Measured_Property> _measured_volumetric_factor;
    std::unique_ptr<Measured_Property> _measured_viscosity;

 public:
    Fluid(){};
    void characterize(int& cells_number);
    void updateProperties(int& term);
    //void calculate(int& term, int& _cellindex);
    const std::string& print() const;
    const int& index() const;
    //sets
    void volumetricFactor(int& _term, int& _cellindex);
    void viscosity(int& _term, int& _cellindex);
    void potential(const int& _term, const int& _cellindex, const double _gravity, const double& depth);

    void density(const int& term, const int& cell_index, const double density){_density[term][cell_index]=density;}
    void pressure(const int& term, const int& cell_index, const double pressure){_pressure[term][cell_index]=pressure;}

    void saturation(const int& term, const int& cell_index, const double saturation){_saturation[term][cell_index]=saturation;}
    
    //gets
    const double& standardConditionsDensity() const {return _standard_conditions_density;};
    const double& pressure             (const int _term, const int _cell_index) const {return _pressure[_term][_cell_index];};             
    const double& density              (const int _term, const int _cell_index) const {return _density[_term][_cell_index];};              
    const double& saturation           (const int _term, const int _cell_index) const {return _saturation[_term][_cell_index];};           
    const double& viscosity            (const int _term, const int _cell_index) const {return _viscosity[_term][_cell_index];};            
    const double& volumetricFactor     (const int _term, const int _cell_index) const {return _volumetric_factor[_term][_cell_index];};    
    const double& potential            (const int _term, const int _cell_index) const {return _potential[_term][_cell_index];};            
    const double& relativePermeability (const int _term, const int _cell_index) const {return _relative_permeability[_term][_cell_index];};

    
};

void Fluid::characterize(int& cells_number){

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();
    std::ostringstream ss = std::ostringstream();
    
    _pressure              = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _density               = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _saturation            = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _viscosity             = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _volumetric_factor     = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _relative_permeability = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));
    _potential             = std::vector<std::vector<double>>(1,std::vector<double>(cells_number));

    myRead(std::string("Please insert the type of fluid "), _type, std::string("Please insert a valid input"));
    
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
    myRead(std::string("Please insert the standard conditions density for the fluid "), _standard_conditions_density, std::string("Please insert a valid input"));
    
    // Initial Conditions for the fluid
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        ss << "Please insert initial pressure for the "<< cellindex+1 << " cell [Pa]";
        myRead(ss.str(), _pressure[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
        
    };
    for(int cellindex=0; cellindex<cells_number; ++cellindex){
        ss << "Please insert initial saturation for the "<< cellindex+1 << " cell [-]";
        myRead(ss.str(), _saturation[0][cellindex], std::string("Please insert a valid input"));
        ss.str("");
        ss.clear();
    };
};

void Fluid::updateProperties(int& term){
    _pressure.push_back             (_pressure[term-1]);
    _potential.push_back            (_potential[term-1]);
    _density.push_back              (_density[term-1]);
    _saturation.push_back           (_saturation[term-1]);
    _viscosity.push_back            (_viscosity[term-1]);
    _volumetric_factor.push_back    (_volumetric_factor[term-1]);
    _relative_permeability.push_back(_relative_permeability[term-1]);
};

void Fluid::volumetricFactor(int& _term, int& _cellindex){
    /*
      Here we need to calculate the restrictions to the flow equations (Volume restriction and Capilarity)

      //Saturation calculation(Volume Restriction)

      //Pressure calculation (Capilarity)

    */
    //Properties Calculation
    _volumetric_factor[_term][_cellindex] =
        _measured_volumetric_factor->interpolate(_pressure[_term][_cellindex]);
};

void Fluid::viscosity(int& _term, int& _cellindex){
    _viscosity[_term][_cellindex] =
        _measured_viscosity->interpolate(_pressure[_term][_cellindex]);
};

const std::string& Fluid::print() const{
    return _type;
};

const int& Fluid::index() const{
    return _index;
};

void Fluid::potential(const int& _term, const int& _cellindex, const double gravity, const double& depth){
    _potential[_term][_cellindex] =
        _pressure[_term][_cellindex] - _density[_term][_cellindex]*gravity*depth;
};

//void {density[_term][_cellindex] = standard_conditions_density;};

#endif /* FLUID_H */
