/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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


#include <aspect/heating_model/latent_heat_lusi.h>

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    LatentHeatLusi<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.n_evaluation_points(),
             ExcMessage ("Heating outputs need to have the same number of entries as the material model inputs."));

      if  (this->introspection().compositional_name_exists("pm_frac"))
        {

          const unsigned int pm_frac_idx = this->introspection().compositional_index_for_name("pm_frac");

          for (unsigned int q=0; q<heating_model_outputs.heating_source_terms.size(); ++q)
            {
              // --- assuming that the entropy_derivative_pressure is 0.0 here (for now)
              //heating_model_outputs.heating_source_terms[q] = 0.0; //material_model_inputs.temperature[q]
                                                                   //* material_model_outputs.densities[q]
                                                          //* material_model_outputs.entropy_derivative_pressure[q]
                                                          //* (material_model_inputs.velocity[q] * material_model_inputs.pressure_gradient[q]);

              // --- pm_frac compo means some presence of asth. partial melt fraction implying an on-going
              //     endothermic reaction hence the plus sign for the entropy_derivative_temperature value
              //     (seems that the matrix assembly code is inconsistent with the energy equation at the
              //     end of section 2.1 of the manual for the LHS term)
              //const double entropy_derivative_temperature=
              //  (material_model_outputs.specific_heat[q]/material_model_inputs.temperature[q])*material_model_inputs.composition[q][pm_frac_idx] ;

              heating_model_outputs.lhs_latent_heat_terms[q] = 0.0; //material_model_outputs.densities[q]
              //                                                 * material_model_inputs.temperature[q] * entropy_derivative_temperature; 
                                                                 //* material_model_outputs.entropy_derivative_temperature[q];

              // --- Apply the latent heat consumption (i.e. effective T cooling ) only if the vertical velo
              //     is upward (implying decompression). The hardcoded 50.0 is the positive
              //     approx. constant value of the delta specific entropy (i.e. S/Kg) for the p.m. of the asth. taken from Connolly 2009
              //     paper (The geodynamic equation of state: What and how, G3 fig. 7 p. ) but it should not
              //     be a constant if we would want to have something more realistic.                
              if (material_model_inputs.velocity[q][dim-1] > 0.0)
                 {

                   // --- Do not use the overloaded vector operator for the dot prod. calculation needlessly
                   //     and assumiing here that the value of
                   //    this->get_gravity_model().gravity_vector(material_model_inputs.position[q])[dim-1]
                   //     is -9.81
                   const double velo_gravity_dot_prod= material_model_inputs.velocity[q][dim-1]
                                                       * this->get_gravity_model().gravity_vector(material_model_inputs.position[q])[dim-1];

                   // ---   ->        ->   ->
                   //      grad P ~= velo . g * density
                   //
                   //      To parametrize the effect of the partial melting of the asth. we simply boost the thermal exp. by 3 orders of
                   //      magnitude in the adabiatic heating-cooling formulat here as a test for a more efficient cooling of the p.m.
                   //      material which we consider as a mixture of the residue (i.e. new OLM material) and p.m. asth pockets to
                   //      decrease its effective T towards the solidus to get higher viscosities (and slower velocities at the same time) .
                   //      With an high (~%15 fraction) p.m. material concentration we then have a more efficient cooling and it is also
                   //      more efficient if the vertical velo value is large (> ~10cm/y)
                   //  
                   heating_model_outputs.heating_source_terms[q] = material_model_outputs.densities[q]
                                                                   * material_model_inputs.temperature[q]
                                                                   * 10e3 * material_model_outputs.thermal_expansion_coefficients[q] 
                                                                   * material_model_inputs.composition[q][pm_frac_idx] * velo_gravity_dot_prod ;
                 }
                 
              // // --- Trying with the equation RHS instead and using dS/dP = -alpha/density relation here but
              // //     since p.m. reaction is endothermic it must cause a cooling of the effective temperature.
              // //     We then omit the minus sign because the p.m. velocity is assumed to be always positive
              // //     upwards and the pressure gradient is considered to be always negative here because of
              // //     the upward displacement of the p.m. compo. The hardcoded 50.0 is the positive
              // //     approx. constant value of the delta entropy for the p.m. of the asth. taken from Connolly 2009
              // //     paper (The geodynamic equation of state: What and how, G3 fig. 7 p. ) but it should not
              // //     be a constant if we would want to have something more realistic.    
              // heating_model_outputs.heating_source_terms[q] = material_model_inputs.composition[q][pm_frac_idx]
              //                                                 * material_model_outputs.thermal_expansion_coefficients[q]
              //                                                 * material_model_inputs.temperature[q]
              //   * (material_model_inputs.velocity[q] * material_model_inputs.pressure_gradient[q]) * 50.0;
 
           }
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(LatentHeatLusi,
                                  "latent heat lusi",
                                  "Implementation of an ad-hoc model specific for the latent heat budget of the Laval U. subduction initiation study")
  }
}
