/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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


#include <aspect/material_model/equation_of_state/multicomponent_incompressible_si.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {

      template <int dim>
      void
      MulticomponentIncompressibleSI<dim>::
      evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
               const unsigned int input_index,
               MaterialModel::EquationOfStateOutputs<dim> &out) const
      {
        // --- Add code to deal with the SI ad-hoc material changes here.
        //    (namely asthenosphere to oceanic lith. mantle and oceanic crust (i.e. basalt) to eclogite)
        std::cout << std::fixed;

        std::cout << "input_index=" << input_index << std::endl;
        //std::cout << "in.temperature[input_index]=" << in.temperature[input_index]<< std::endl;
        //std::cout << "in.pressure[input_index]=" << in.pressure[input_index]<< std::endl;

        //std::vector<double> composition (in.composition[input_index]);

        //for (unsigned int c=0; c<in.composition[input_index].size(); ++c)
        //{
        //    std::cout << "\nbef modif: this->introspection().name_for_compositional_index(c)="
        //              << this->introspection().name_for_compositional_index(c) << std::endl;
        //    std::cout << "bef modif: in.composition[input_index][c]=" << in.composition[input_index][c] << std::endl;
        //}

        const unsigned int asth_mtl_idx= this->introspection().
           compositional_index_for_name(asthenospheric_mantle_nid);

        const unsigned int oc_lith_mtl_idx= this->introspection().
           compositional_index_for_name(oceanic_lithospheric_mantle_nid);

        //const unsigned int oc_crust_basalts_idx= this->introspection().
        //   compositional_index_for_name(oceanic_crust_basalts_nid);

        // --- VERY UGLY TRICK done here to bypass the const MaterialModel::MaterialModelInputs<dim>
        //     to be able to modify the composition on-the-fly
        std::vector<double>* composition_ref=
            (std::vector<double, std::allocator<double> >*)&in.composition[input_index];

        //std::cout << "composition_ref=" <<
        //printf("in.composition=%xi\n",in.composition);
        //printf("composition_ref=%xi\n",composition_ref);

        // --- local copy of oc. ltih. mantle compo.
        //const double oc_lith_mtl_compo_tmp= (*composition_ref)[oc_lith_mtl_idx];

        (*composition_ref)[asth_mtl_idx]= 0.5;
        (*composition_ref)[oc_lith_mtl_idx]= 0.25;

        // --- NOTE: The asth. material now becomes oc. lith. mantle material at this input_index
        //           BUT we need to keep the value of the oc. lith. mantle (could be
        //           > 0.0) that was already defined and add the asth. compo value to
        //           to it
        (*composition_ref)[oc_lith_mtl_idx] += in.composition[input_index][asth_mtl_idx]; //(*composition_ref)[asth_mtl_idx]

        // --- Avoid negative compo value
        (*composition_ref)[oc_lith_mtl_idx]= std::max(0.0,(*composition_ref)[oc_lith_mtl_idx]);

        // --- Can now set the asth. compo to 0.0 since this value has been
        //     completely transferred to the oc. lith. mantle compo value.
        //     (mass conservation principle)
        (*composition_ref)[asth_mtl_idx]= 0.0; //oc_lith_mtl_compo_tmp;

        //for (unsigned int c=0; c<in.composition[input_index].size(); ++c)
        //{
        //    std::cout << "\naft modif: this->introspection().name_for_compositional_index(c)="
        //              << this->introspection().name_for_compositional_index(c) << std::endl;
        //    std::cout << "aft modif: in.composition[input_index][c]=" << in.composition[input_index][c] << std::endl;
        //}

        std::exit(0);

        // --- now resume with the usual MulticomponentIncompressible super class code exec.
        MulticomponentIncompressible<dim>::evaluate((const MaterialModel::MaterialModelInputs<dim>)in,input_index,out);
      }


      template <int dim>
      bool
      MulticomponentIncompressibleSI<dim>::
      is_compressible () const
      {
        return MulticomponentIncompressible<dim>::is_compressible();
      }



      template <int dim>
      void
      MulticomponentIncompressibleSI<dim>::declare_parameters (ParameterHandler &prm,
                                                             const double default_thermal_expansion)
      {
         MulticomponentIncompressible<dim>::declare_parameters(prm,default_thermal_expansion);
      }


      template <int dim>
      void
      MulticomponentIncompressibleSI<dim>::parse_parameters (ParameterHandler &prm,
                                                           const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {
        MulticomponentIncompressible<dim>::parse_parameters(prm,expected_n_phases_per_composition);
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
#define INSTANTIATE(dim) \
  template class MulticomponentIncompressibleSI<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}
