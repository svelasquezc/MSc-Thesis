#ifndef FLUID_H
#define FLUID_H

#include <vector>
#include <memory>
#include <string>
#include <iostream>

#include "Measured_Property.h"

class Fluid{
 private:
    int index;
    std::string type;
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
    Fluid(){};
    void Characterize(){};
    
};

#endif /* FLUID_H */
