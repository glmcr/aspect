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
      //     (i.e. rock type transformation depending on
      //     the dynamic and thermodynamic conditions)
      const unsigned int asth_mtl_idx=
	this->introspection().compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID);

      const unsigned int oc_lith_mtl_idx=
	this->introspection().compositional_index_for_name(LITHOSPHERIC_MANTLE_NID);

      const unsigned int oc_crust_idx=
	this->introspection().compositional_index_for_name(OCEANIC_CRUST_NID);

      const unsigned int olm_asth_hybrid_idx=
	this->introspection().compositional_index_for_name(OLM_ASTH_HYBRID_NID);

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

               //// --- Determine if the material is locally undergoing decompression or
	       ////     compression at this ith evaluation point. We are using the sign
	       ////     of the scalar product of the velocity and the pressure_gradient
	       ////     for doing so in an almost self-consistent manner. We have decompression
	       ////     if this scalar product is negative and compression if it is positive
	       ////     regardless of the respective orientations of the velocity and pressure
	       ////     gradient (i.e. we do not assume that their vertical components are
	       ////     always larger than their respective horizontal components)
	       //const bool decompression= (in.velocity[i] * in.pressure_gradient[i] < 0.0);

               // --- 1st ad-hoc material change (rock type transformation):
               //     asthenosphere becomes basaltic oceanic crust, lithospheric
	       //     mantle or an hybrid between the two depending on the temperature,
	       //     pressure, velocity direction and pressure gradient direction.
               if (in.temperature[i] <= LAB_TEMPERATURE_IN_KELVINS)
                 {

                   // --- Determine if the material is locally undergoing decompression or
                   //     compression at this ith evaluation point. We are using the sign
                   //     of the scalar product of the velocity and the pressure_gradient
                   //     for doing so in an almost self-consistent manner. We have decompression
                   //     if this scalar product is negative and compression if it is positive
                   //     regardless of the respective orientations of the velocity and pressure
                   //     gradient (i.e. we do not assume that their vertical components are
                   //     always larger than their respective horizontal components)
                   //const bool decompression= (in.velocity[i] * in.pressure_gradient[i] < 0.0);

                   // --- Get the asthenospheric composition value (if any, could be 0.0)
                   const double ast_reaction_term= in.composition[i][asth_mtl_idx];

                   // --- The asth. is not allowed to exist at T < LAB_TEMPERATURE_IN_KELVINS
		   //     It will become 0.0 for the next time step.
                   out.reaction_terms[i][asth_mtl_idx]= -ast_reaction_term;

                   // --- Also we possibly need to transform the hybrid material (if any) that was created
                   //     with the asthenospheric mantle in compression conditions to the "normal"
                   //     lithospheric mantle OR oceanic crust
                   const double hyb_reaction_term= in.composition[i][olm_asth_hybrid_idx];

		   // --- Within the parametrized ad-hoc SI context here the decompression implies that
		   //     the asthenospheric mantle is undergoing partial fusion (even if the material is cooling)
		   //     and the basaltic crust is then created (i.e. appears suddenly) if the local pressure is
		   //     <= MOHO_PRESSURE_IN_PASCALS or the lithospheric mantle (harzburgite) is created here
		   //     otherwise if the local pressure is > MOHO_PRESSURE_IN_PASCALS.
                   if (in.pressure[i] <= MOHO_PRESSURE_IN_PASCALS)
                     {

                         // --- Upwelling asth. mantle and hybrid material transform to oc. crust via
                         //     the out.reaction_terms
                         out.reaction_terms[i][oc_crust_idx]=
                           ast_reaction_term + hyb_reaction_term;

                         //// --- hybrid material transformed to lithospheric mantle
                         //out.reaction_terms[i][oc_lith_mtl_idx]= hyb_2_lmt_reaction_term;

                     }
		   else if (in.pressure[i] <= OLM_MAX_PRESSURE_IN_PASCALS)
                     {

                      // --- Here the SI ad-hoc parametrization implies that the asth. transform to
                      //     the lithos. mantle (harzburgite) via the out.reaction_terms data vector.
	              //     (In reality, the harzburgite is the solid residue of the partial
	              //      fusion of the asthenospheric mantle which is dynamically accreted
	              //      to the solid sides of the oceanic ridge or supra-subduction zone
		      //      oceanic lithosphere). We also need to add the hybrid material
                      //      since it is supposed to behave like asthenospheric material here.
                      out.reaction_terms[i][oc_lith_mtl_idx]=
                        ast_reaction_term + hyb_reaction_term;
                     }
                    else
                      {
                        // ---  pressure > OLM_MAX_PRESSURE_IN_PASCALS.
                        //     No partial fusion, asth. becomes the hybrid material here.
                        out.reaction_terms[i][olm_asth_hybrid_idx]= ast_reaction_term;

                        // --- Since there is no loss of the hybrid material here
                        //     then its reaction term need to be zero.
                        hyb_reaction_term= 0.0;

                      } // --- end inner if-else block

                      // --- Need to remove the hybrid material reaction term from
                      //     itself here because it was possibly used in one of the
                      //     first two above transformations (to oc. crust OR to lith. mantle).
                      //     But hyb_reaction_term could be zero if the 3rd transformation
                      //     need to be done.
                      out.reaction_terms[i][olm_asth_hybrid_idx]= -hyb_reaction_term

		   //else // --- Here the material is undergoing cooling and compression -> no partial fusion
		   //  {
		   //     // --- The asthenospheric mantle here becomes an hybrid rock material having the thermal cond.
		   //     //     of the lithospheric mantle and having the asthenospheric mantle viscous flow law. We
                   //     //     also suppose that this hybrid material has the same geochemistry as the
                   //     //     asth. mantle (Note: this decrease of th thermal cond. seems to be going against the
                   //     //     fact that the thermal cond. should increase when a material cools, not decrease. Here
                   //     //     this ad-hoc parametrization is applied to mimic that the convective asth. mantle becomes
                   //     //     non-convective when its temperature becomes < 1573K (1300C).
                   //     out.reaction_terms[i][olm_asth_hybrid_idx]= ast_reaction_term;
                   //
		   //  } // --- end outer if-else block

                   // --- And finally the asthenospheric mantle composition (concentration) need to
                   //     to become zero at this evaluation point i at the next time step because its
                   //     concentration has now been transferred to another rock material (either oc. crust
                   //     OR lith. mantle OR hybrid material) at the same evaluation point (Note the usage
                   //     of the minus sign here on the rhs for the assignation of the out.reaction_terms[i][asth_mtl_idx]
                   //     on the lhs).
                   //out.reaction_terms[i][asth_mtl_idx]= -ast_reaction_term;

	       }  // --- end if (in.temperature[i] <= LAB_TEMPERATURE_IN_KELVINS) block

               //else // --- here T > LAB_TEMPERATURE_IN_KELVINS
               //  {
               //
               //    const double lmt_2_ast_reaction_term= in.composition[i][oc_lith_mtl_idx];
               //
               //    // --- Now the lithospheric mantle transform to asthenosphere at evaluation point i
               //    //     (but only if in.composition[i][oc_lith_mtl_idx] > 0.0)
               //    out.reaction_terms[i][asth_mtl_idx]= lmt_2_ast_reaction_term;
               //       //in.composition[i][oc_lith_mtl_idx]; // * inv_current_time_step;
               //
               //    // --- And the lithospheric mantle disappear (if it was not 0.0!!)
               //    //     at evaluation point i
               //    out.reaction_terms[i][oc_lith_mtl_idx]= -lmt_2_ast_reaction_term;
               //    //   -in.composition[i][oc_lith_mtl_idx]; // * inv_current_time_step;

               //  } // --- inner if-else block

               // --- TODO: implement the lithospheric mantle -> basalt (one way only) material change
               //           when T>= 873 (moho T)

               // --- 2nd and more complicated ad-hoc material change:
               //     basaltic material of the oceanic crust becomes eclogite.

	  } // --- inner for loop block on  in.n_evaluation_points()
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

