#ifndef INTERFLUID_INTERACTION_H
#define INTERFLUID_INTERACTION_H

#include "Fluid.h"

class Interfluid_Interaction : protected Value_Reader {
 private:
    static int _count_of_interfluid_interactions;
    int _index;
    std::weak_ptr<Fluid> _reference_fluid;
    std::weak_ptr<Fluid> _wetting_fluid;
    std::weak_ptr<Fluid> _non_wetting_fluid;
    std::unique_ptr<Measured_Property> _measured_reference_relative_permeability;
    std::unique_ptr<Measured_Property> _measured_principal_relative_permeability;
    std::unique_ptr<Measured_Property> _measured_capillary_pressure;

 public:

    Interfluid_Interaction(){};
 Interfluid_Interaction(const int index) : _index(index) {};
    
    const double irreducibleSaturation() const
    {
        return _measured_principal_relative_permeability->begin()->first;
    };
    
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

    const std::shared_ptr<Fluid> referenceFluid()  const {return _reference_fluid.lock();}
    const std::shared_ptr<Fluid> wettingFluid()    const {return _wetting_fluid.lock();}
    const std::shared_ptr<Fluid> nonWettingFluid() const {return _non_wetting_fluid.lock();}

    void add(std::vector<std::shared_ptr<Fluid>>& characterized_fluids){

	int counter;
	int reference;
	int wetting_fliud;
	int non_wetting_fliud;

        std::ostringstream ref_name = std::ostringstream();
        std::ostringstream ref_value = std::ostringstream();

	if(characterized_fluids.size() >= 2){

	    std::cout << "Please select the Reference Fluid: " << std::endl;
	    for (counter=0; counter<characterized_fluids.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), reference, std::string("Please insert a valid index"));
                if(reference>0 && reference<=characterized_fluids.size()){
                    _reference_fluid = characterized_fluids[reference-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            std::cout << "Please select the Wetting Fluid: " << std::endl;
	    for (counter=0; counter<characterized_fluids.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), wetting_fliud, std::string("Please insert a valid index"));
                if(wetting_fliud>0 && wetting_fliud<=characterized_fluids.size()){
                    _wetting_fluid = characterized_fluids[wetting_fliud-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            std::cout << "Please select the Non Wetting Fluid: " << std::endl;
	    for (counter=0; counter<characterized_fluids.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_fluids[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), non_wetting_fliud, std::string("Please insert a valid index"));
                if(non_wetting_fliud>0 && non_wetting_fliud<=characterized_fluids.size()){
                    _non_wetting_fluid = characterized_fluids[non_wetting_fliud-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            ref_name << _reference_fluid.lock()->print() << " Saturation";

            ref_value << _reference_fluid.lock()->print() << "Reference Relative Permeability";
        
            _measured_reference_relative_permeability =
                std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

            ref_value.str("");
            ref_value.clear();

            for (auto fluid : characterized_fluids){
                if(fluid->principal()){
                    ref_value << fluid->print() << " Relative Permeability to " << _reference_fluid.lock()->print();

                    _measured_principal_relative_permeability =
                        std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

                    ref_value.str("");
                    ref_value.clear();

                    ref_value << _non_wetting_fluid.lock()->print() << "-" << _wetting_fluid.lock()->print() << "Capillary Pressure";

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

    void addFromFile(std::ifstream& interaction_reader, std::vector<std::shared_ptr<Fluid>>& characterized_fluids){

        std::string element;
        int counter;
	int reference;
	int wetting_fluid;
	int non_wetting_fluid;

        std::ostringstream ref_name = std::ostringstream();
        std::ostringstream ref_value = std::ostringstream();

        while(interaction_reader >> element){

            if(element == "REFERENCE_FLUID"){

                interaction_reader >> reference;

                if(reference < 1 && reference>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Reference fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                }else{
                    _reference_fluid = characterized_fluids[reference-1];

                    ref_name << _reference_fluid.lock()->print() << " Saturation";

                    ref_value << _reference_fluid.lock()->print() << "Reference Relative Permeability";
        
                    _measured_reference_relative_permeability =
                        std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

                    ref_value.str("");
                    ref_value.clear();

                    for (auto fluid : characterized_fluids){
                        if(fluid->principal()){
                            ref_value << fluid->print() << " Relative Permeability to " << _reference_fluid.lock()->print();

                            _measured_principal_relative_permeability =
                                std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));

                            ref_value.str("");
                            ref_value.clear();

                        };

                    };
                        
                };

            }else if(element == "WETTING_FLUID"){
                
                interaction_reader >> wetting_fluid;

                if(wetting_fluid < 1 && wetting_fluid>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Wetting_Fluid fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                }else{
                    _wetting_fluid = characterized_fluids[wetting_fluid-1];
                };
                
            }else if(element == "NON_WETTING_FLUID"){
                interaction_reader >> non_wetting_fluid;

                if(non_wetting_fluid < 1 && non_wetting_fluid>characterized_fluids.size()){
                    std::stringstream fluid_index_err = std::stringstream();
                    fluid_index_err << "Non_Wetting_Fluid fluid index must be between 1 and "<<characterized_fluids.size()<<std::endl;
                    throw std::out_of_range(fluid_index_err.str());
                    
                }else{
                    _non_wetting_fluid = characterized_fluids[non_wetting_fluid-1];
                };
            }else if(element == "REFERENCE_RELATIVE_PERMEABILITY"){
                
                _measured_reference_relative_permeability->readFromFile(interaction_reader);
                
            }else if(element == "PRINCIPAL_RELATIVE_PERMEABILITY"){
                
                _measured_principal_relative_permeability->readFromFile(interaction_reader);
                
            }else if(element == "CAPILLARY_PRESSURE"){
                    
                ref_value << _non_wetting_fluid.lock()->print() << "-" << _wetting_fluid.lock()->print() << "Capillary Pressure";

                _measured_capillary_pressure =
                    std::make_unique<Measured_Property>(Measured_Property(ref_name.str(),ref_value.str()));
                
                _measured_capillary_pressure->readFromFile(interaction_reader);

                break;
            };
        };
    };
};
#endif /* INTERFLUID_INTERACTION_H */
