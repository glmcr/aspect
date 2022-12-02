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

#include <aspect/material_model/visco_plastic_si.h>
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
    bool
    ViscoPlasticSI<dim>::
    is_yielding(const MaterialModelInputs<dim> &in) const
    {
      return ViscoPlastic<dim>::is_yielding(in) ;//plastic_yielding;
    }


    template <int dim>
    void
    ViscoPlasticSI<dim>::
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

               // --- 1st ad-hoc material change (rock type transformation):
               //     asthenosphere becomes lithospheric mantle
               //     if the T of the former at location i is at or under LAB_TEMPERATURE_IN_KELVINS
               if (in.temperature[i] <= LAB_TEMPERATURE_IN_KELVINS)
                 {

                   // ---
                   const double ast_2_lmt_reaction_term= in.composition[i][asth_mtl_idx];

                   if (in.pressure[i] <= MOHO_PRESSURE_IN_PASCALS)
                     {

                       // --- asth. transform to oc. crust via the out.reaction_terms
                       out.reaction_terms[i][oc_crust_idx]= ast_2_lmt_reaction_term;

                       // --- And the asthenosphere composition (concentration) will become zero at
                       //     evaluation point i at the next time step because of the usage of
                       //     this reaction term (note the minus sign here)
                       //out.reaction_terms[i][asth_mtl_idx]= -ast_2_lmt_reaction_term;

                     }
                   else
                     {

                       // --- asth. transform to Lithos. mantle via the out.reaction_terms
                       out.reaction_terms[i][oc_lith_mtl_idx]= ast_2_lmt_reaction_term;

                       // --- And the asthenosphere composition (concentration) will become zero at
                       //     evaluation point i at the next time step because of the usage of
                       //     this reaction term (note the minus sign here)
                       //out.reaction_terms[i][asth_mtl_idx]= -ast_2_lmt_reaction_term;
                       //   -in.composition[i][asth_mtl_idx]; //* inv_current_time_step;
                     }

                     // --- And the asthenosphere composition (concentration) will become zero at
                     //     evaluation point i at the next time step because of the usage of
                     //     this reaction term (note the minus sign here)
                     out.reaction_terms[i][asth_mtl_idx]= -ast_2_lmt_reaction_term;
                 }
               else // --- here T > LAB_TEMPERATURE_IN_KELVINS
                 {

                   const double lmt_2_ast_reaction_term= in.composition[i][oc_lith_mtl_idx];

                   // --- Now the lithospheric mantle transform to asthenosphere at evaluation point i
                   //     (but only if in.composition[i][oc_lith_mtl_idx] > 0.0)
                   out.reaction_terms[i][asth_mtl_idx]= lmt_2_ast_reaction_term;
                      //in.composition[i][oc_lith_mtl_idx]; // * inv_current_time_step;

                   // --- And the lithospheric mantle disappear (if it was not 0.0!!)
                   //     at evaluation point i
                   out.reaction_terms[i][oc_lith_mtl_idx]= -lmt_2_ast_reaction_term;
                   //   -in.composition[i][oc_lith_mtl_idx]; // * inv_current_time_step;

                 } // --- inner if-else block

               // --- TODO: implement the lithospheric mantle -> basalt (one way only) material change
               //           when T>= 873 (moho T)

               // --- 2nd and more complicated ad-hoc material change:
               //     basaltic material of the oceanic crust becomes eclogite.

	  } // --- for loop block
        } // --- outer if block
      } // --- method block


    template <int dim>
    bool
    ViscoPlasticSI<dim>::
    is_compressible () const
    {
      return ViscoPlastic<dim>::is_compressible(); //equation_of_state.is_compressible();
    }



    template <int dim>
    double ViscoPlasticSI<dim>::
    get_min_strain_rate () const
    {
      return ViscoPlastic<dim>::get_min_strain_rate(); //rheology->min_strain_rate;
    }



    template <int dim>
    void
    ViscoPlasticSI<dim>::declare_parameters (ParameterHandler &prm)
    {
       ViscoPlastic<dim>::declare_parameters(prm);
    }



    template <int dim>
    void
    ViscoPlasticSI<dim>::parse_parameters (ParameterHandler &prm)
    {

       //this->get_pcout() << std::endl <<"ViscoPlasticSI::parse_parameters(), Hello world!!" << std::endl;

       ViscoPlastic<dim>::parse_parameters(prm);

       //this->get_pcout() << std::endl <<"ViscoPlasticSI::parse_parameters(), done with ViscoPlastic<dim>::parse_parameters(prm)" << std::endl;
    }



    template <int dim>
    void
    ViscoPlasticSI<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      ViscoPlastic<dim>::create_additional_named_outputs(out);
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(ViscoPlasticSI,
                                   "visco plastic si",
                                   "Temporary implementation of a class that inherits from the ViscoPlastic class "
                                   "and which allows to add ad-hoc material changes using the reaction rates. This is "
                                   "done in relation to a Ph. D. study on subduction initiation.")
  }
}

