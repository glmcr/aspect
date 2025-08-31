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

      // --- Need to use the normal generic evaluate method of the ViscoPlastic super class
      ViscoPlastic<dim>::evaluate(in,out);

      // // --- apply the LUSI ad-hoc code that deals with
      // //     the variation of the thermal expansivity and diffusivity
      // //     of the SI materials.
      //ViscoPlastic<dim>::apply_lusi_thermal_stuff(in,out);
      
      // // --- ***IMPORTANT*** Need to add one to the indices
      // //     since the material models properties arrays are always using
      // //     a background field at index 0 so the compositions indices
      // //     are shifted by one in the material models properties arrays
      // //const unsigned int asth_mtl_idx= this->introspection().
      // //                   compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID) + 1;

      // //const unsigned int oc_lith_mtl_idx= this->introspection().
      // //                   compositional_index_for_name(LITHOSPHERIC_MANTLE_NID) + 1;

      // //const unsigned int oc_crust_idx= this->introspection().
      // //                   compositional_index_for_name(OCEANIC_CRUST_NID) + 1;
      
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

