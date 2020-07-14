#ifndef INTERPHASE_INTERACTION_H
#define INTERPHASE_INTERACTION_H

#include "Phase.h"

class Interphase_Interaction : protected Value_Reader {
private:
    static int _count_of_interphase_interactions;
    int _index;
    std::weak_ptr<Phase> _reference_phase;
    std::weak_ptr<Phase> _wetting_phase;
    std::weak_ptr<Phase> _non_wetting_phase;
    std::unique_ptr<Property_Measure> _reference_relative_permeability_measure;
    std::unique_ptr<Property_Measure> _main_relative_permeability_measure;
    std::unique_ptr<Property_Measure> _capillary_pressure_measure;

public:

    Interphase_Interaction(){};
    Interphase_Interaction(const int index) : _index(index) {};
    
    const double irreducibleSaturation() const
    {
        return _main_relative_permeability_measure->begin()->first;
    };

    const double maximumMainRelativePermeability() const
    {
        return (--_main_relative_permeability_measure->end())->second;
    };
    
    double referenceRelativePermeability(const double saturation) const
    {
	return _reference_relative_permeability_measure->interpolate(saturation);
    };

    double mainRelativePermeability(const double saturation) const
    {
	return _main_relative_permeability_measure->interpolate(saturation);
    };

    double capillaryPressure(const double saturation) const
    {
	return _capillary_pressure_measure->interpolate(saturation);
    };

    const std::shared_ptr<Phase> referencePhase()  const {return _reference_phase.lock();}
    const std::shared_ptr<Phase> wettingPhase()    const {return _wetting_phase.lock();}
    const std::shared_ptr<Phase> nonWettingPhase() const {return _non_wetting_phase.lock();}

    void add(std::vector<std::shared_ptr<Phase>>& characterized_phases){

	int counter;
	int reference;
	int wetting_phase;
	int non_wetting_phase;

        std::ostringstream ref_name = std::ostringstream();
        std::ostringstream ref_value = std::ostringstream();

	if(characterized_phases.size() >= 2){

	    std::cout << "Please select the Reference Phase: " << std::endl;
	    for (counter=0; counter<characterized_phases.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), reference, std::string("Please insert a valid index"));
                if(reference>0 && reference<=characterized_phases.size()){
                    _reference_phase = characterized_phases[reference-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            std::cout << "Please select the Wetting Phase: " << std::endl;
	    for (counter=0; counter<characterized_phases.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), wetting_phase, std::string("Please insert a valid index"));
                if(wetting_phase>0 && wetting_phase<=characterized_phases.size()){
                    _wetting_phase = characterized_phases[wetting_phase-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            std::cout << "Please select the Non Wetting Phase: " << std::endl;
	    for (counter=0; counter<characterized_phases.size(); ++counter){
		std::cout << (counter+1) << ". " << characterized_phases[counter]->print() << std::endl;
	    };
	    
            while(true){
                myRead(std::string(""), non_wetting_phase, std::string("Please insert a valid index"));
                if(non_wetting_phase>0 && non_wetting_phase<=characterized_phases.size()){
                    _non_wetting_phase = characterized_phases[non_wetting_phase-1];
                    break;
                }else{
                    std::cout << "Please insert an index inside the range" << std::endl;
                }
            
            };

            ref_name << _reference_phase.lock()->print() << " Saturation";

            ref_value << _reference_phase.lock()->print() << "Reference Relative Permeability";
        
            _reference_relative_permeability_measure =
                std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));

            ref_value.str("");
            ref_value.clear();

            for (auto phase : characterized_phases){
                if(phase->main()){
                    ref_value << phase->print() << " Relative Permeability to " << _reference_phase.lock()->print();

                    _main_relative_permeability_measure =
                        std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));

                    ref_value.str("");
                    ref_value.clear();

                    ref_value << _non_wetting_phase.lock()->print() << "-" << _wetting_phase.lock()->print() << "Capillary Pressure";

                    _capillary_pressure_measure =
                        std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
                };

            };

            _reference_relative_permeability_measure->readMe();
            _main_relative_permeability_measure->readMe();
            _capillary_pressure_measure->readMe();

	}else{
	    
	    std::cout << "It is not possible to add an Interphase interaction with only one phase characterized."
                      << std::endl;
	    
	};
    };

    void addFromFile(std::ifstream& interaction_reader, std::vector<std::shared_ptr<Phase>>& characterized_phases){

        std::string element;
        int counter;
	int reference;
	int wetting_phase;
	int non_wetting_phase;

        std::ostringstream ref_name = std::ostringstream();
        std::ostringstream ref_value = std::ostringstream();

        while(interaction_reader >> element){

            if(element == "REFERENCE_PHASE"){

                interaction_reader >> reference;

                if(reference < 1 && reference>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Reference phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                }else{
                    _reference_phase = characterized_phases[reference-1];

                    ref_name << _reference_phase.lock()->print() << " Saturation";

                    ref_value << _reference_phase.lock()->print() << "Reference Relative Permeability";
        
                    _reference_relative_permeability_measure =
                        std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));

                    ref_value.str("");
                    ref_value.clear();

                    for (auto phase : characterized_phases){
                        if(phase->main()){
                            ref_value << phase->print() << " Relative Permeability to " << _reference_phase.lock()->print();

                            _main_relative_permeability_measure =
                                std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));

                            ref_value.str("");
                            ref_value.clear();

                        };

                    };
                        
                };

            }else if(element == "WETTING_PHASE"){
                
                interaction_reader >> wetting_phase;

                if(wetting_phase < 1 && wetting_phase>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Wetting_Phase phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                }else{
                    _wetting_phase = characterized_phases[wetting_phase-1];
                };
                
            }else if(element == "NON_WETTING_PHASE"){
                interaction_reader >> non_wetting_phase;

                if(non_wetting_phase < 1 && non_wetting_phase>characterized_phases.size()){
                    std::stringstream phase_index_err = std::stringstream();
                    phase_index_err << "Non_Wetting_Phase phase index must be between 1 and "<<characterized_phases.size()<<std::endl;
                    throw std::out_of_range(phase_index_err.str());
                    
                }else{
                    _non_wetting_phase = characterized_phases[non_wetting_phase-1];
                };
            }else if(element == "REFERENCE_RELATIVE_PERMEABILITY"){
                
                _reference_relative_permeability_measure->readFromFile(interaction_reader);
                
            }else if(element == "MAIN_RELATIVE_PERMEABILITY"){
                
                _main_relative_permeability_measure->readFromFile(interaction_reader);
                
            }else if(element == "CAPILLARY_PRESSURE"){
                    
                ref_value << _non_wetting_phase.lock()->print() << "-" << _wetting_phase.lock()->print() << "Capillary Pressure";

                _capillary_pressure_measure =
                    std::make_unique<Property_Measure>(Property_Measure(ref_name.str(),ref_value.str()));
                
                _capillary_pressure_measure->readFromFile(interaction_reader);

                break;
            };
        };
    };
};
#endif /* INTERPHASE_INTERACTION_H */
