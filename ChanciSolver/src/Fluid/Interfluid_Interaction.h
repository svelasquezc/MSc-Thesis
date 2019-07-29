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

    void add(const int fluids_quantity, std::vector<std::shared_ptr<Fluid>>& MyFluids){

	int counter;
	int reference;
	int wetting_fliud;
	int non_wetting_fliud;

    std::ostringstream ref_name = std::ostringstream();
    std::ostringstream ref_value = std::ostringstream();

	if(fluids_quantity >= 2){

	    std::cout << "Please select the Reference Fluid: " << std::endl;
	    for (counter=0; counter<MyFluids.size(); ++counter){
		std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
	    };
	    
        while(true){
            myRead(std::string(""), reference, std::string("Please insert a valid index"));
            if(reference>0 && reference<=MyFluids.size()){
                _reference_fluid = MyFluids[reference-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };

        std::cout << "Please select the Wetting Fluid: " << std::endl;
	    for (counter=0; counter<MyFluids.size(); ++counter){
		std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
	    };
	    
        while(true){
            myRead(std::string(""), wetting_fliud, std::string("Please insert a valid index"));
            if(wetting_fliud>0 && wetting_fliud<=MyFluids.size()){
                _wetting_fluid = MyFluids[wetting_fliud-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };

        std::cout << "Please select the Non Wetting Fluid: " << std::endl;
	    for (counter=0; counter<MyFluids.size(); ++counter){
		std::cout << (counter+1) << ". " << MyFluids[counter]->print() << std::endl;
	    };
	    
        while(true){
            myRead(std::string(""), non_wetting_fliud, std::string("Please insert a valid index"));
            if(non_wetting_fliud>0 && non_wetting_fliud<=MyFluids.size()){
                _non_wetting_fluid = MyFluids[non_wetting_fliud-1];
                break;
            }else{
                std::cout << "Please insert an index inside the range" << std::endl;
            }
            
        };

        ref_name << _reference_fluid->print() << " Saturation"

        ref_value << _reference_fluid->print() << " Relative Permeability";
        
        _measured_reference_relative_permeability =
            std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

        ref_value.str("");
        ref_value.clear();

        for (auto fluid : MyFluids){
            if(fluid->principal()){
                ref_value << fluid->print() << " Relative Permeability to " << _reference_fluid->print();

                _measured_principal_relative_permeability =
                    std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

                ref_value.str("");
                ref_value.clear();

                ref_value << fluid->print() << "-" << _reference_fluid->print() << "Capillary Pressure";

                _measured_capillary_pressure =
                    std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
            };

        };

        _measured_reference_relative_permeability->readMe();
        _measured_principal_relative_permeability->readMe();
        _measured_capillary_pressure->readMe();

	}else{
	    
	    std::cout << "It is not possible to add an Interfluid interaction with only one fluid characterized."
                  << std::endl;
	    
	};
    };
};

#endif /* INTERFLUID_INTERACTION_H */
