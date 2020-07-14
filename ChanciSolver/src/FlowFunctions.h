#ifndef FLOW_FUNCTIONS_H
#define FLOW_FUNCTIONS_H

#include "Initial_Conditions.h"
#include "Database.h"

using namespace Database;

void calculateInteractions(const int& term, const int& cell_index, Phase& phase){

    double capillary_pressure=0;
    
    for(auto& interphase_interaction : added_interphase_interactions )
        {
            auto reference_phase = interphase_interaction->referencePhase();
            auto wetting_phase = interphase_interaction->wettingPhase();
            auto non_wetting_phase = interphase_interaction->nonWettingPhase();
	    
            if(phase.index() == reference_phase->index()){
		
                phase.relativePermeability(term, cell_index,
                                           interphase_interaction->referenceRelativePermeability(phase.saturation(term, cell_index)));
                capillary_pressure = interphase_interaction->capillaryPressure(phase.saturation(term, cell_index));
	    
                if(phase.index() == wetting_phase->index()){
                    phase.pressure(term, cell_index,
                                   non_wetting_phase->pressure(term, cell_index) - capillary_pressure);
                }else{
                    phase.pressure(term, cell_index,
                                   wetting_phase->pressure(term, cell_index) + capillary_pressure);
                };
            };
        };
};

double calculateBaker(const int& term, const int& cell_index){

    double accumulated_saturation = 0;
    double accumulated_main_relative_permeability=0;
    double irreducible_saturation = 0;
    double mobile_saturation=0;
    double interpolated_main_relative_permeability=0;
    double maximum_main_relative_permeability=0;
    

    for(auto& interphase_interaction : added_interphase_interactions )
        {
            auto reference_phase = interphase_interaction->referencePhase();
	    
            irreducible_saturation = interphase_interaction->irreducibleSaturation();
	    
            interpolated_main_relative_permeability = interphase_interaction->
                mainRelativePermeability(reference_phase->saturation(term, cell_index));

            mobile_saturation = reference_phase->saturation(term, cell_index) - irreducible_saturation;
	    
            accumulated_saturation += mobile_saturation;

            accumulated_main_relative_permeability += mobile_saturation*interpolated_main_relative_permeability;

            maximum_main_relative_permeability = std::max(maximum_main_relative_permeability, interphase_interaction->maximumMainRelativePermeability());
            
        };
    
    if(accumulated_saturation == 0){
        return maximum_main_relative_permeability;
    }else{
        return accumulated_main_relative_permeability/accumulated_saturation;
    };
    
};

void calculateProperties(const int& term, const std::shared_ptr<Cell>& cell, Rock& rock){

    const auto cell_index = cell->index();

    double remaining_saturation = 1.0;

    for(auto phase : characterized_phases){

        if(phase->main()){
            rock.porosity(term, cell_index, phase->pressure(term, cell_index));
            //rock.porosity(term, cell_index, phase->pressure(term, cell_index));
        }else{
            remaining_saturation = remaining_saturation - phase->saturation(term, cell_index);

            calculateInteractions(term, cell_index, *phase);
        };

    };

    for(auto phase : characterized_phases){

        if(phase->main()){
            phase->saturation(term, cell_index, remaining_saturation);
            if(Initial_Conditions::phases_quantity>=2){
                phase->relativePermeability(term, cell_index, calculateBaker(term, cell_index));
            }else{
                phase->relativePermeability(term, cell_index, 1.0);
            };
        };

        phase->formationVolumeFactor(term, cell_index, phase->pressure(term, cell_index));
        phase->viscosity(term, cell_index, phase->pressure(term, cell_index));
        //Calculate density

        double density_contribution = phase->standardConditionsDensity();
        
        for(auto& equilibrium_relationship : added_equilibrium_relationships ){
            
            if(equilibrium_relationship->contributingPhase()->index() == phase->index()){
                
                const auto receiving = equilibrium_relationship->receivingPhase();
                // I suppose this should interpolate the partition coefficient
                // to the contributing phase pressure
                equilibrium_relationship->partitionCoefficient(term,cell_index);

                density_contribution = density_contribution + 
                    equilibrium_relationship->partitionCoefficient(term,cell_index) *
                    receiving->standardConditionsDensity();
                
            };
        };

        density_contribution = density_contribution/phase->formationVolumeFactor(term, cell_index);

        phase->density(term, cell_index, density_contribution);
        
        phase->potential(term, cell_index, gravity, cell->depth());
    };        
};

double calculateAccumulation(const int& term, Phase& phase, const std::shared_ptr<Cell>& cell, Rock& rock){

    double past_contribution=0;
    double current_contribution=0;

    //std::string N="N";
    //std::string K="K";

    const int cell_index = cell->index();
    
    for(auto& equilibrium_relationship : added_equilibrium_relationships ){
        
        if(equilibrium_relationship->receivingPhase()->index() == phase.index()){
            
            const auto contributing = equilibrium_relationship->contributingPhase();

            double past_coef = equilibrium_relationship->partitionCoefficient(term-1,cell_index);
            
            past_contribution = past_contribution +
                past_coef * (rock.porosity(term-1,cell_index) * contributing->saturation(term-1,cell_index)
                             / contributing->formationVolumeFactor(term-1,cell_index));

            double curr_coef = equilibrium_relationship->partitionCoefficient(term,cell_index);
            
            current_contribution = current_contribution +
                curr_coef * (rock.porosity(term,cell_index) * contributing->saturation(term,cell_index)
                             / contributing->formationVolumeFactor(term,cell_index));
            
        };
    };

    double accumulation = (cell->volume()/Initial_Conditions::timedelta) *
        (((rock.porosity(term,cell_index)*phase.saturation(term,cell_index)
           /phase.formationVolumeFactor(term,cell_index)) + current_contribution)-
         ((rock.porosity(term-1,cell_index)*phase.saturation(term-1,cell_index)
           /phase.formationVolumeFactor(term-1,cell_index)) + past_contribution));
    
    return accumulation;
};

double calculateFlow(const int& term, Phase& phase, const Mesh& mesh, const std::shared_ptr<Cell>& cell, const std::shared_ptr<Face>& face, Rock& rock){
    
    auto harmonicAverage = [](double cell_property, double neighbor_property){
                                                                              return 1.0/((1.0/cell_property) + (1.0/neighbor_property));
    };

    double flow=0;
    
    int direction = face->orientation();

    auto neighbor_cell = face->neighborCell().lock();
    
    const auto neighbor_index = neighbor_cell->index();
    const auto cell_index = cell->index();

    const auto neighbor_axis = neighbor_cell->numeration3D()[direction];
    const auto axis = cell->numeration3D()[direction];

    double length_delta = mesh.thickness(direction,axis);
    double neighbor_length_delta = mesh.thickness(direction,neighbor_axis);

    double porous_volume = rock.porosity(term, cell_index)*cell->volume();
    double neighbor_porous_volume = rock.porosity(term, neighbor_axis)*neighbor_cell->volume();

    // here is corrected the shape_factor calculation, review the preconceptual schema...
    double shape_factor = face->area()*rock.absolutePermeability(term, cell_index)[direction]/length_delta;
    double neighbor_shape_factor = face->area()*rock.absolutePermeability(term, neighbor_index)[direction]/neighbor_length_delta;

    double face_volumetric_factor = (phase.formationVolumeFactor(term, cell_index)*porous_volume +
                                     phase.formationVolumeFactor(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_viscosity = (phase.viscosity(term, cell_index)*porous_volume +
                             phase.viscosity(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_shape_factor = 2.0*harmonicAverage(shape_factor, neighbor_shape_factor);

    double face_relative_permeability=1;
    
    //Upwind!!!!
    if(phase.potential(term, cell_index) >= phase.potential(term, neighbor_index)){
        face_relative_permeability = phase.relativePermeability(term, cell_index);
    }else{
        face_relative_permeability = phase.relativePermeability(term, neighbor_index);
    };

    double face_transmissivity = face_shape_factor * face_relative_permeability /
        (face_volumetric_factor * face_viscosity);

    flow = flow + face_transmissivity*(phase.potential(term, neighbor_index) - phase.potential(term, cell_index));

    for(auto& equilibrium_relationship : added_equilibrium_relationships ){
        
        if(equilibrium_relationship->receivingPhase()->index() == phase.index()){

            const auto contributing = equilibrium_relationship->contributingPhase();
            
            double face_volumetric_factor =
                (contributing->formationVolumeFactor(term, cell_index)*
                 porous_volume +
                 contributing->formationVolumeFactor(term, neighbor_index)*
                 neighbor_porous_volume) /
                (porous_volume + neighbor_porous_volume);

            double face_viscosity = (contributing->viscosity(term, cell_index)*porous_volume +
                                     contributing->viscosity(term, neighbor_index)*neighbor_porous_volume) /
                (porous_volume + neighbor_porous_volume);

            //Upwind!!!!
            if(contributing->potential(term, cell_index) >= contributing->potential(term, neighbor_index)){
                face_relative_permeability = contributing->relativePermeability(term, cell_index);
            }else{
                face_relative_permeability = contributing->relativePermeability(term, neighbor_index);
            };

            double face_partition_coefficient =
                (equilibrium_relationship->partitionCoefficient(term, cell_index)*porous_volume +
                 equilibrium_relationship->partitionCoefficient(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);
            
            double face_transmissivity = face_shape_factor * face_relative_permeability /
                (face_volumetric_factor * face_viscosity);

            flow = flow + face_partition_coefficient*face_transmissivity*
                (contributing->potential(term, neighbor_index) - contributing->potential(term, cell_index));
        };
    };

    neighbor_cell.reset();
    
    return flow;
};

#endif /*FLOW_FUNCTIONS_H*/
