#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <memory>
#include <string>
#include <sstream>

#include "Measured_Property.h"

class Fluid : protected Value_Reader{
    
 private:
    
    int index;
    std::string type;
    std::vector<std::vector<double>> pressure;
    std::vector<std::vector<double>> density;
    std::vector<std::vector<double>> saturation;
    std::vector<std::vector<double>> viscosity;
    std::vector<std::vector<double>> volumetric_factor;
    std::vector<std::vector<double>> potential;
    std::vector<std::vector<double>> relative_permeability;
    double standard_conditions_density;
    std::unique_ptr<Measured_Property> measured_volumetric_factor;
    std::unique_ptr<Measured_Property> measured_viscosity;

 public:
    Fluid(){
        // Creation of secure pointers for properties
        
        measured_volumetric_factor = std::make_unique<Measured_Property>();
        measured_viscosity = std::make_unique<Measured_Property>();
    };
    void characterize(int& cells_number);
    void updateProperties(int& term);
    //void calculate(int& term, int& _cellindex);
    void calculateVolumetricFactor(int& _term, int& _cellindex);
    void calculateViscosity(int& _term, int& _cellindex);
    const std::string print();
    int getIndex();
};

void Fluid::characterize(int& _cells_number){

    std::stringstream ss = std::stringstream();
    
    pressure              = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    density               = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    saturation            = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    viscosity             = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    volumetric_factor     = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    relative_permeability = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));
    potential             = std::vector<std::vector<double>>(1,std::vector<double>(_cells_number));

    myRead(std::string("Please insert the type of fluid\n"), type, std::string("Please insert a valid input"));
    
    // Reading of measured properties (PVT Table)
    measured_volumetric_factor->readMe();
    measured_viscosity->readMe();
    
    // Reading of standard conditions density
    myRead(std::string("Please insert the standard conditions density for the fluid\n"), standard_conditions_density, std::string("Please insert a valid input"));
    
    // Initial Conditions for the fluid
    for(int cellindex=0; cellindex<_cells_number; ++cellindex){
        ss << "Please insert initial pressure for the "<< cellindex+1 << "cell" <<std::endl;
        myRead(ss.str(), pressure[0][cellindex], std::string("Please insert a valid input"));
        ss.flush();
        
        ss << "Please insert initial saturation for the "<< cellindex+1 << "cell" <<std::endl;
        myRead(ss.str(), saturation[0][cellindex], std::string("Please insert a valid input"));
        ss.flush();
    };
};

void Fluid::updateProperties(int& term){
    pressure.push_back(pressure[term-1]);
    density.push_back(density[term-1]);
    saturation.push_back(saturation[term-1]);
    viscosity.push_back(viscosity[term-1]);
    volumetric_factor.push_back(volumetric_factor[term-1]);
    relative_permeability.push_back(relative_permeability[term-1]);
    potential.push_back(potential[term-1]);
};

void Fluid::calculateVolumetricFactor(int& _term, int& _cellindex){
    /*
      Here we need to calculate the restrictions to the flow equations (Volume restriction and Capilarity)

      //Saturation calculation(Volume Restriction)

      //Pressure calculation (Capilarity)

    */
    //Properties Calculation
    volumetric_factor[_term][_cellindex] =
        measured_volumetric_factor->interpolate(pressure[_term][_cellindex]);
};

void Fluid::calculateViscosity(int& _term, int& _cellindex){
    viscosity[_term][_cellindex] =
        measured_viscosity->interpolate(pressure[_term][_cellindex]);
};

const std::string Fluid::print(){
    return type;
};

int Fluid::getIndex(){
    return index;
}

//void {density[_term][_cellindex] = standard_conditions_density;};

#endif /* FLUID_H */
