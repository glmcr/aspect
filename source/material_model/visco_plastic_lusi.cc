/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/

#include <aspect/material_model/visco_plastic_lusi.h>
#include <aspect/utilities.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/newton.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>
//#include <aspect/lusi_composition.h>

namespace aspect
{
  namespace MaterialModel
  {

    template <int dim>
    void
    ViscoPlasticLUSI<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {

      ViscoPlastic<dim>::evaluate(in,out);

      // --- ***IMPORTANT*** Need to add one to the indices
      //     since the material models properties arrays are always using
      //     a background field at index 0 so the compositions indices
      //     are shifted by one in the material models properties arrays
      //const unsigned int asth_mtl_idx= this->introspection().
      //                   compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID) + 1;

      //const unsigned int oc_lith_mtl_idx= this->introspection().
      //                   compositional_index_for_name(LITHOSPHERIC_MANTLE_NID) + 1;

      //const unsigned int oc_crust_idx= this->introspection().
      //                   compositional_index_for_name(OCEANIC_CRUST_NID) + 1;

      const unsigned int ssz_oc_crust_idx= this->introspection().
	compositional_index_for_name(LUSIComposition<dim>::SSZ_OCEANIC_CRUST_NID) + 1;
      //compositional_index_for_name(Particle::Property::LUSIComposition<dim>::SSZ_OCEANIC_CRUST_NID) + 1;
      
      // --- Only apply the ad-hoc material changes if the simulator initialization
      //     is done.
      if  (this->simulator_is_past_initialization() && this->get_timestep_number() > 0 )
        {
          const ComponentMask volumetric_compositions = this->rheology->get_volumetric_composition_mask();
	  
          const EquationOfState::MulticomponentIncompressible<dim> eos_cref= this->equation_of_state_constref;

          const double reference_T_cref= eos_cref.reference_T_constref;
          
          const std::vector<double> densities_cref= eos_cref.densities_constref;
          const std::vector<double> thermal_expansivities_cref= eos_cref.thermal_expansivities_constref;

          //EquationOfStateOutputs<dim> eos_out();
	  
          std::vector<double> densities_local(eos_cref.densities_constref);
          std::vector<double> thermal_expansivities_local(eos_cref.thermal_expansivities_constref);
          
          const double reference_temperature = reference_T_cref;

          // --- Loop through all requested points
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {

               //this->get_pcout() << std::endl << "ViscoPlasticLUSI::execute: in.temperature[i]= " << in.temperature[i] << std::endl ;

	       const std::vector<double> volume_fractions= MaterialUtilities::
		     compute_composition_fractions(in.composition[i], volumetric_compositions);

               if (in.temperature[i] > THERMAL_EXP_LOW_T_IN_K_THRESHOLD && in.temperature[i] < THERMAL_EXP_UPP_T_IN_K_THRESHOLD)
                 {

                   //this->get_pcout() << std::endl << "ViscoPlasticLUSI::execute: in.temperature[i]= "
                   //                  << in.temperature[i] << std::endl ;

                   // If adiabatic heating is used, the reference temperature used to calculate density should be the adiabatic
                   // temperature at the current position. This definition is consistent with the Extended Boussinesq Approximation.
                   //const double reference_temperature = (this->include_adiabatic_heating()
                   //                                      ?
                   //                                      this->get_adiabatic_conditions().temperature(in.position[i])
                   //                                      :
                   //                                      reference_T_cref);

                   //this->get_pcout() << "ViscoPlasticLUSI:: reference_temperature=" << reference_temperature << std::endl ;

                   const double thExpFact= 1.0 +
                       (in.temperature[i]-THERMAL_EXP_LOW_T_IN_K_THRESHOLD)*THERMAL_EXP_T_IN_K_THRD_FACT;

                   // thermal_expansivities_local[asth_mtl_idx]=
                   //      thExpFact*thermal_expansivities_cref[asth_mtl_idx];

		   //this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_cref[asth_mtl_idx]= "
                   //                  << thermal_expansivities_cref[asth_mtl_idx] << std::endl;
		   //this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_local[asth_mtl_idx]= "
                   //                  << thermal_expansivities_local[asth_mtl_idx] << std::endl;

                   // densities_local[asth_mtl_idx]= densities_cref[asth_mtl_idx] *
                   //    (1.0 - thermal_expansivities_local[asth_mtl_idx] * (in.temperature[i] - reference_temperature));

                   // thermal_expansivities_local[oc_lith_mtl_idx]=
                   //     thExpFact*thermal_expansivities_cref[oc_lith_mtl_idx];

		   //this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_cref[oc_lith_mtl_idx]= "
                   //                  << thermal_expansivities_cref[oc_lith_mtl_idx] << std::endl;
		   //this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_local[oc_lith_mtl_idx]= "
                   //                  << thermal_expansivities_local[oc_lith_mtl_idx] << std::endl << std::endl;

                   // densities_local[oc_lith_mtl_idx]= densities_cref[oc_lith_mtl_idx] *
                   //  (1.0 - thermal_expansivities_local[oc_lith_mtl_idx]* (in.temperature[i] - reference_temperature));

		   // --- Update thermal expansivitires and densities accordinglym (for all compos)
                   for (unsigned int cmp=0; cmp < volume_fractions.size(); ++cmp)
		   {
		      thermal_expansivities_local[cmp]= thExpFact * thermal_expansivities_cref[cmp];

                      densities_local[cmp]= densities_cref[cmp] *
                          (1.0 - thermal_expansivities_local[cmp] * (in.temperature[i] - reference_temperature));		      
		   }

                   out.densities[i]=
		      MaterialUtilities::average_value (volume_fractions, densities_local, MaterialUtilities::arithmetic);

                   out.thermal_expansion_coefficients[i]=
                      MaterialUtilities::average_value (volume_fractions, thermal_expansivities_local, MaterialUtilities::arithmetic);

                   //if (in.temperature[i] > THERMAL_EXP_LOW_T_IN_K_THRESHOLD && in.temperature[i] < 1600) {
                   //  this->get_pcout() << std::endl << "ViscoPlasticLUSI::execute: in.temperature[i]= "
                   //                  << in.temperature[i] << std::endl ;
                   //  this->get_pcout() << "asth_mtl_idx=" << asth_mtl_idx << std::endl ;
                   //  this->get_pcout() << "oc_lith_mtl_idx=" << oc_lith_mtl_idx << std::endl ;
                   //  this->get_pcout() << "oc_crust_idx=" << oc_crust_idx << std::endl ;
                   //  this->get_pcout() << "ViscoPlasticLUSI:: reference_temperature="
                   //                     << reference_temperature << std::endl ;
                   //  this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_cref[asth_mtl_idx]= "
                   //                  << thermal_expansivities_cref[asth_mtl_idx] << std::endl;
                   //  this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_local[asth_mtl_idx]= "
                   //                  << thermal_expansivities_local[asth_mtl_idx] << std::endl;
                   //  this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_cref[oc_lith_mtl_idx]= "
                   //                    << thermal_expansivities_cref[oc_lith_mtl_idx] << std::endl;
                   //  this->get_pcout() << "ViscoPlasticLUSI:: thermal_expansivities_local[oc_lith_mtl_idx]= "
                   //                  << thermal_expansivities_local[oc_lith_mtl_idx] << std::endl << std::endl;
                   //}
		 } // --- if block for thermal exp. dependance on T

	       // --- Now update the thermal cond. via the thermal_diffusivities.
	       if ( in.temperature[i] < THERMAL_EXP_UPP_T_IN_K_THRESHOLD) {

		 double thermal_diffusivity= 0.0;

		 // --- Local copy of the this->thermal_diffusivities vector
		 //     to be able to (artificially) increase the thermal_diffusivities
		 //     of the material(s) that could be subjected to hydrothermal
		 //     convection or circulation (mainly the newly created oc. crust at
		 //     an oceanic ridge). This is an ad-hoc brute-force parametrization.
                 std::vector<double> thermal_diffusivities_local(this->thermal_diffusivities);

                 // --- Increase the ssz oc. crust thermal diff using the strain rate as a
		 //     proxy for the hydrothermal convection intensity near the ridge. A
		 //     large strain rate implies a vigourous hydrothermal convection.
		 //     (it is at the limit between viscous and plastic strains domains)
                 //     It helps to localize the heat dissipation at the ridge which itself
		 //     promotes the strain localization (or vice-versa?) with cracking-faulting.
  
                 // --- Get the strain rate that was used for the viscosity calculation at i
                 const double min_strain_rate= this->get_min_strain_rate();
		 
                 const double edot_ii= 
		     std::max(std::sqrt(std::max(-second_invariant(deviator(in.strain_rate[i])), 0.)),min_strain_rate);

		 // --- Apply the increasing factor using the strain rate and the min strain rate values

		 // --- Need to use an hardcoded 3e7 factor for the min_strain_rate
		 //     in the denom part of the exponential. Also limit the strain rate
		 //     max to MAX_STRAIN_RATE_LOCAL to avoid producing a factor > 30.  
                 thermal_diffusivities_local[ssz_oc_crust_idx] *= exp( ( std::min(edot_ii,MAX_STRAIN_RATE_LOCAL) - min_strain_rate)/(min_strain_rate*3e7) );
		 
		 for (unsigned int cmp=0; cmp < volume_fractions.size(); ++cmp)
		 {
		   thermal_diffusivity += volume_fractions[cmp] * thermal_diffusivities_local[cmp];
		   //thermal_diffusivity += volume_fractions[cmp] * this->thermal_diffusivities[cmp];
                 }

		 // --- NOTE: We assume here that the reference T is 273K
		 //     Limit the thDiffFactor between 1.0 and 0.45
		 const double thDiffFactor=
		   std::min(1.0,std::max(1.0 - THERMAL_DIFF_T_IN_K_FACT*(in.temperature[i] - reference_temperature), 0.45));
		 
		 thermal_diffusivity *= thDiffFactor ;

		 out.thermal_conductivities[i]= thermal_diffusivity * out.specific_heat[i] * out.densities[i];
		 
	       } // --- end if block for thermal cond. dependance on T
	       
            } // --- inner for loop
        } // --- outer if block
    } // --- method block
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlasticLUSI,
                                   "visco plastic lusi",
                                   "Temporary implementation of a class that inherits from the ViscoPlastic class "
                                   "and which allows to do ad-hoc materical model modifications related to G. Mercier Ph.D. thesis at U. Laval")
  }
}

