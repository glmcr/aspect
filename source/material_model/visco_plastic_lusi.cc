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

      // --- Now take care of the ad-hoc material changes
      //     (i.e. rock type transformation depending only
      //     on T or on both p and T)
      const unsigned int asth_mtl_idx= this->introspection().
           compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID);

      const unsigned int oc_lith_mtl_idx= this->introspection().
           compositional_index_for_name(LITHOSPHERIC_MANTLE_NID);

      const unsigned int oc_crust_idx= this->introspection().
           compositional_index_for_name(OCEANIC_CRUST_NID);

      // --- Only apply the ad-hoc material changes if the simulator initialization
      //     is done.
      if  (this->simulator_is_past_initialization() && this->get_timestep_number() > 0 )
        {

          // --- No need to use this->get_timestep() here
          //     (It's only relevant for the reaction rates it seems)
          //const double inv_current_time_step= 1.0/this->get_timestep();

          // --- Loop through all requested points
          for (unsigned int i=0; i < in.n_evaluation_points(); ++i)
            {

               //this->get_pcout() << "ViscoPlasticSI::execute: in.temperature[i]= " << in.temperature[i] << std::endl;
               if (in.temperature[i] > THERMAL_EXP_LOW_T_IN_K_THRESHOLD && in.temperature[i] < THERMAL_EXP_UPP_T_IN_K_THRESHOLD)
                 {

		   
		 }
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

