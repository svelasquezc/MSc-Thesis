#ifndef INTERFLUID_INTERACTION_H
#define INTERFLUID_INTERACTION_H

#include "Fluid.h"

class Interfluid_Interaction : protected Value_Reader {
 private:
    static int _count_of_interfluid_interactions;
    int _index;
    std::shared_ptr<Fluid> _reference_fluid;
    std::shared_ptr<Fluid> _wetting_fluid;
    std::shared_ptr<Fluid> _non_wetting_fluid;
    std::unique_ptr<Measured_Property> _measured_reference_relative_permeability;
    std::unique_ptr<Measured_Property> _measured_principal_relative_permeability;
    std::unique_ptr<Measured_Property> _measured_capillary_pressure;

 public:

    double referenceRelativePermeability(const double saturation) const
    {
	return _measured_reference_relative_permeability->interpolate(saturation);
    };

    double principalRelativePermeability(const double saturation) const
    {
	return _measured_principal_relative_permeability->interpolate(saturation);
    };

    double capillaryPressure(const double saturation) const
    {
	return _measured_capillary_pressure->interpolate(saturation);
    };

    const std::shared_ptr<Fluid>& referenceFluid()  const {return _reference_fluid;}
    const std::shared_ptr<Fluid>& wettingFluid()    const {return _wetting_fluid;}
    const std::shared_ptr<Fluid>& nonWettingFluid() const {return _non_wetting_fluid;}
};

#endif /* INTERFLUID_INTERACTION_H */
