#ifndef FLOW_FUNCTIONS_H
#define FLOW_FUNCTIONS_H

#include "Initial_Conditions.h"
#include "Database.h"

using namespace Database;

void calculateInteractions(const int& term, const int& cell_index, Fluid& fluid){

    double capillary_pressure=0;
    
    for(auto& interfluid_interaction : added_interfluid_interactions )
        {
            auto reference_fluid = interfluid_interaction->referenceFluid();
            auto wetting_fluid = interfluid_interaction->wettingFluid();
            auto non_wetting_fluid = interfluid_interaction->nonWettingFluid();
	    
            if(fluid.index() == reference_fluid->index()){
		
                fluid.relativePermeability(term, cell_index,
                                           interfluid_interaction->referenceRelativePermeability(fluid.saturation(term, cell_index)));
                capillary_pressure = interfluid_interaction->capillaryPressure(fluid.saturation(term, cell_index));
	    
                if(fluid.index() == wetting_fluid->index()){
                    fluid.pressure(term, cell_index,
                                   non_wetting_fluid->pressure(term, cell_index) - capillary_pressure);
                }else{
                    fluid.pressure(term, cell_index,
                                   wetting_fluid->pressure(term, cell_index) + capillary_pressure);
                };
            };
        };
};

double calculateBaker(const int& term, const int& cell_index){

    double accumulated_saturation = 0;
    double accumulated_principal_relative_permeability=0;
    double irreducible_saturation = 0;
    double mobile_saturation=0;
    double interpolated_principal_relative_permeability=0;
    double maximum_principal_relative_permeability=0;
    

    for(auto& interfluid_interaction : added_interfluid_interactions )
        {
            auto reference_fluid = interfluid_interaction->referenceFluid();
	    
            irreducible_saturation = interfluid_interaction->irreducibleSaturation();
	    
            interpolated_principal_relative_permeability = interfluid_interaction->
                principalRelativePermeability(reference_fluid->saturation(term, cell_index));

            mobile_saturation = reference_fluid->saturation(term, cell_index) - irreducible_saturation;
	    
            accumulated_saturation += mobile_saturation;

            accumulated_principal_relative_permeability += mobile_saturation*interpolated_principal_relative_permeability;

            maximum_principal_relative_permeability = std::max(maximum_principal_relative_permeability, interfluid_interaction->maximumPrincipalRelativePermeability());
            
        };
    
    if(accumulated_saturation == 0){
        return maximum_principal_relative_permeability;
    }else{
        return accumulated_principal_relative_permeability/accumulated_saturation;
    };
    
};

void calculateProperties(const int& term, const std::shared_ptr<Cell>& cell, Rock& rock){

    const auto cell_index = cell->index();

    double remaining_saturation = 1.0;

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            rock.porosity(term, cell_index, fluid->pressure(term, cell_index));
            //rock.porosity(term, cell_index, fluid->pressure(term, cell_index));
        }else{
            remaining_saturation = remaining_saturation - fluid->saturation(term, cell_index);

            calculateInteractions(term, cell_index, *fluid);
        };

    };

    for(auto fluid : characterized_fluids){

        if(fluid->principal()){
            fluid->saturation(term, cell_index, remaining_saturation);
            if(Initial_Conditions::fluids_quantity>=2){
                fluid->relativePermeability(term, cell_index, calculateBaker(term, cell_index));
            }else{
                fluid->relativePermeability(term, cell_index, 1.0);
            };
        };

        fluid->volumetricFactor(term, cell_index, fluid->pressure(term, cell_index));
        fluid->viscosity(term, cell_index, fluid->pressure(term, cell_index));
        //Calculate density

        double density_contribution = fluid->standardConditionsDensity();
        
        for(auto& equilibrium_relation : added_equilibrium_relations ){
            
            if(equilibrium_relation->contributorFluid()->index() == fluid->index()){
                
                const auto receiver = equilibrium_relation->receiverFluid();
                // I suppose this should interpolate the partition coefficient
                // to the contributor fluid pressure
                equilibrium_relation->partitionCoefficient(term,cell_index);

                density_contribution = density_contribution + 
                    equilibrium_relation->partitionCoefficient(term,cell_index) *
                    receiver->standardConditionsDensity();
                
            };
        };

        density_contribution = density_contribution/fluid->volumetricFactor(term, cell_index);

        fluid->density(term, cell_index, density_contribution);
        
        fluid->potential(term, cell_index, gravity, cell->depth());
    };        
};

double calculateAccumulation(const int& term, Fluid& fluid, const std::shared_ptr<Cell>& cell, Rock& rock){

    double past_contribution=0;
    double current_contribution=0;

    //std::string N="N";
    //std::string K="K";

    const int cell_index = cell->index();
    
    for(auto& equilibrium_relation : added_equilibrium_relations ){
        
        if(equilibrium_relation->receiverFluid()->index() == fluid.index()){
            
            const auto contributor = equilibrium_relation->contributorFluid();

            double past_coef = equilibrium_relation->partitionCoefficient(term-1,cell_index);
            
            past_contribution = past_contribution +
                past_coef * (rock.porosity(term-1,cell_index) * contributor->saturation(term-1,cell_index)
                             / contributor->volumetricFactor(term-1,cell_index));

            double curr_coef = equilibrium_relation->partitionCoefficient(term,cell_index);
            
            current_contribution = current_contribution +
                curr_coef * (rock.porosity(term,cell_index) * contributor->saturation(term,cell_index)
                             / contributor->volumetricFactor(term,cell_index));
            
        };
    };

    double accumulation = (cell->volume()/Initial_Conditions::timedelta) *
        (((rock.porosity(term,cell_index)*fluid.saturation(term,cell_index)
           /fluid.volumetricFactor(term,cell_index)) + current_contribution)-
         ((rock.porosity(term-1,cell_index)*fluid.saturation(term-1,cell_index)
           /fluid.volumetricFactor(term-1,cell_index)) + past_contribution));
    
    return accumulation;
};

double calculateFlow(const int& term, Fluid& fluid, const Mesh& mesh, const std::shared_ptr<Cell>& cell, const std::shared_ptr<Face>& face, Rock& rock){
    
    auto harmonicAverage = [](double cell_property, double neighbor_property){
        return 1.0/((1.0/cell_property) + (1.0/neighbor_property));
    };

    double flow=0;
    
    int direction = face->orientation();

    auto neighbor_cell = face->neighbor().lock();
    
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

    double face_volumetric_factor = (fluid.volumetricFactor(term, cell_index)*porous_volume +
                                     fluid.volumetricFactor(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_viscosity = (fluid.viscosity(term, cell_index)*porous_volume +
                             fluid.viscosity(term, neighbor_index)*neighbor_porous_volume) /
        (porous_volume + neighbor_porous_volume);

    double face_shape_factor = 2.0*harmonicAverage(shape_factor, neighbor_shape_factor);

    double face_relative_permeability=1;
    
    //Upwind!!!!
    if(fluid.potential(term, cell_index) >= fluid.potential(term, neighbor_index)){
        face_relative_permeability = fluid.relativePermeability(term, cell_index);
    }else{
        face_relative_permeability = fluid.relativePermeability(term, neighbor_index);
    };

    double face_transmissivity = face_shape_factor * face_relative_permeability /
        (face_volumetric_factor * face_viscosity);

    flow = flow + face_transmissivity*(fluid.potential(term, neighbor_index) - fluid.potential(term, cell_index));

    for(auto& equilibrium_relation : added_equilibrium_relations ){
        
        if(equilibrium_relation->receiverFluid()->index() == fluid.index()){

            const auto contributor = equilibrium_relation->contributorFluid();
            
            double face_volumetric_factor = (contributor->volumetricFactor(term, cell_index)*porous_volume +
                                             contributor->volumetricFactor(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);

            double face_viscosity = (contributor->viscosity(term, cell_index)*porous_volume +
                                     contributor->viscosity(term, neighbor_index)*neighbor_porous_volume) /
                (porous_volume + neighbor_porous_volume);

            //Upwind!!!!
            if(contributor->potential(term, cell_index) >= contributor->potential(term, neighbor_index)){
                face_relative_permeability = contributor->relativePermeability(term, cell_index);
            }else{
                face_relative_permeability = contributor->relativePermeability(term, neighbor_index);
            };

            double face_partition_coefficient =
                (equilibrium_relation->partitionCoefficient(term, cell_index)*porous_volume +
                 equilibrium_relation->partitionCoefficient(term, neighbor_index)*neighbor_porous_volume) / (porous_volume + neighbor_porous_volume);
            
            double face_transmissivity = face_shape_factor * face_relative_permeability /
                (face_volumetric_factor * face_viscosity);

            flow = flow + face_partition_coefficient*face_transmissivity*
                (contributor->potential(term, neighbor_index) - contributor->potential(term, cell_index));
        };
    };

    neighbor_cell.reset();
    
    return flow;
};

#endif /*FLOW_FUNCTIONS_H*/
