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

#ifndef _aspect_material_model_visco_plastic_lusi_h
#define _aspect_material_model_visco_plastic_lusi_h

#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/visco_plastic.h>
#include<deal.II/fe/component_mask.h>

#include <aspect/particle/property/lusi_composition.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    static constexpr const double THERMAL_EXP_UPP_T_IN_K_THRESHOLD= Particle::Property::ASTH_OLM_HYBRID_MAT_TEMP_THESHOLD_KELVINS;

    template <int dim>
    class ViscoPlasticLUSI : public MaterialModel::ViscoPlastic<dim>
    {
      public:

      void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                    MaterialModel::MaterialModelOutputs<dim> &out) const override;
        // /**
        //  * Return whether the model is compressible or not.  Incompressibility
        //  * does not necessarily imply that the density is constant; rather, it
        //  * may still depend on temperature or pressure. In the current
        //  * context, compressibility means whether we should solve the continuity
        //  * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
        //  * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
        //  *
        //  * This material model is incompressible.
        //  */
        // bool is_compressible () const override;

        // static
        // void
        // declare_parameters (ParameterHandler &prm);

        // void
        // parse_parameters (ParameterHandler &prm) override;

        // void
        // create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        // double get_min_strain_rate() const;

        // /**
        //  * A function that returns whether the material is plastically yielding at
        //  * the given pressure, temperature, composition, and strain rate.
        //  *
        //  * @deprecated: Use the other function with this name instead, which allows
        //  * to pass in more general input variables.
        //  */
        // //DEAL_II_DEPRECATED
        // //bool
        // //is_yielding (const double pressure,
        // //             const double temperature,
        // //             const std::vector<double> &composition,
        // //             const SymmetricTensor<2,dim> &strain_rate) const;

        // /**
        //  * A function that returns whether the material is plastically
        //  * yielding at the given input variables (pressure, temperature,
        //  * composition, strain rate, and so on).
        //  */
        // bool
        // is_yielding (const MaterialModelInputs<dim> &in) const;

        // // ---
        // //inline
        static constexpr const char* ASTHENOSPHERIC_MANTLE_NID= "asthenosphere";

        // // ---
        // //inline
        static constexpr const char* LITHOSPHERIC_MANTLE_NID= "oceanicLithMantle";

        // // ---
        static constexpr const char* OCEANIC_CRUST_NID= "oceanicCrust";

        static constexpr const double THERMAL_EXP_LOW_T_IN_K_THRESHOLD= 500.0;
      //static constexpr const double THERMAL_EXP_UPP_T_IN_K_THRESHOLD= 2000.0;
      
      // use the ASTH_OLM_HYBRID_MAT_TEMP_THESHOLD_KELVINS T to define the upper limit
      // of the application of the dependance of the thermal exp. on T. This is to
      // keep the asthenosphere at the constant 52.5 thermal cond as in Butler 2017
      //  static constexpr const double THERMAL_EXP_UPP_T_IN_K_THRESHOLD= Particle::Property::ASTH_OLM_HYBRID_MAT_TEMP_THESHOLD_KELVINS;
      
        static constexpr const double THERMAL_EXP_T_IN_K_THRD_FACT=
          0.22/(THERMAL_EXP_UPP_T_IN_K_THRESHOLD-THERMAL_EXP_LOW_T_IN_K_THRESHOLD);
      
        // --- Factor to apply to the thermal diff. to parametrize its dependance on
        //     the temperature. Note that the thermal cond. is already depending on
        //     the material density (which itself depends on the temperature) for the
        //     visco-plastic material but we want to have the thermal cond. depending
        //     on both the density AND the thermal diff.
        //     NOTE: We assume here that the reference T is 273K
        //           1573K is the T threshold for the asth. vs lith. mantle limit.
        static constexpr const double THERMAL_DIFF_T_IN_K_FACT= 1.0/(2.0*(1573.0-273.0));
      
      //private:
      //
      //  /**
      //   * Pointer to the object used to compute the rheological properties.
      //   * In this case, the rheology in question is visco(elasto)plastic. The
      //    * object contains functions for parameter declaration and parsing,
      //   * and further functions that calculate viscosity and viscosity
      //   * derivatives. It also contains functions that create and fill
      //   * additional material model outputs, specifically plastic outputs.
      //   * The rheology itself is a composite rheology, and so the object
      //   * contains further objects and/or pointers to objects that provide
      //   * functions and parameters for all subordinate rheologies.
      //   */
      //  std::unique_ptr<Rheology::ViscoPlastic<dim>> rheology;
      //
      //  std::vector<double> thermal_diffusivities;
      //
      //  /**
      //   * Whether to use user-defined thermal conductivities instead of thermal diffusivities.
      //   */
      //  bool define_conductivities;
      //
      //  std::vector<double> thermal_conductivities;
      //
      //  /**
      //   * Object for computing the equation of state.
      //   */
      //EquationOfState::MulticomponentIncompressibleLUSI<dim> equation_of_state;
      //
      //  /**
      //   * Object that handles phase transitions.
      //   */
      //  MaterialUtilities::PhaseFunction<dim> phase_function;

    };
  }
}

#endif
