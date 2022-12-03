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

#ifndef _aspect_material_model_visco_plastic_si_h
#define _aspect_material_model_visco_plastic_si_h

#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/visco_plastic.h>
//#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
//#include <aspect/material_model/rheology/visco_plastic.h>
#include <aspect/material_model/utilities_si.h>

#include<deal.II/fe/component_mask.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    template <int dim>
    class ViscoPlasticSI : public MaterialModel::ViscoPlastic<dim> //, public ::aspect::SimulatorAccess<dim>
    {
      public:

        void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                      MaterialModel::MaterialModelOutputs<dim> &out) const override;

        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         *
         * This material model is incompressible.
         */
        bool is_compressible () const override;

        static
        void
        declare_parameters (ParameterHandler &prm);

        void
        parse_parameters (ParameterHandler &prm) override;

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

        double get_min_strain_rate() const;

        /**
         * A function that returns whether the material is plastically yielding at
         * the given pressure, temperature, composition, and strain rate.
         *
         * @deprecated: Use the other function with this name instead, which allows
         * to pass in more general input variables.
         */
        //DEAL_II_DEPRECATED
        //bool
        //is_yielding (const double pressure,
        //             const double temperature,
        //             const std::vector<double> &composition,
        //             const SymmetricTensor<2,dim> &strain_rate) const;

        /**
         * A function that returns whether the material is plastically
         * yielding at the given input variables (pressure, temperature,
         * composition, strain rate, and so on).
         */
        bool
        is_yielding (const MaterialModelInputs<dim> &in) const;

        // ---
        //inline
        static constexpr const char* ASTHENOSPHERIC_MANTLE_NID= "asthenosphere";

        // ---
        //inline
        static constexpr const char* LITHOSPHERIC_MANTLE_NID= "oceanicLithMantle";

        // --- Oceanic lithos.  <-> asthenoshperic hybrid material
        //     (asth. flow law but OLM thermal cond. to parametrize
        //     the conversion of convecting asth. mantle to lithospheric
        //     mantle by cooling only i.e. not by partial fusion)
        //     It does not exists at the beginning and is produced
        //     when the asthenospheric material cools because of
        //     the downward migration of the 1573K (1300C) isotherm
        //     and the adiabatic heating is positive (implying
        //     very slow downward material advection which would
        //     be slower than the isotherm downward migration itself)
        static constexpr const char* OLM_ASTH_HYBRID_NID= "olmAsthHybrid";
      
        // ---
        static constexpr const char* OCEANIC_CRUST_NID= "oceanicCrust";

        // --- Lithosphere <-> asthenosphere T boundary
        static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1573.0;

        // ---
        static constexpr const double GPA_2_PASCALS= 1e9;

        static constexpr const double KBARS_2_GPA= 0.1;

        static constexpr const double KBARS_2_PASCALS= KBARS_2_GPA*GPA_2_PASCALS;

        // --- ~2.5Kb is the lithostatic pressure at the oceanic moho depth (~7km)
        //     It is used to separate the transformation of the upwelling asthenosphere
        //     to lithospheric mantle (p>2.5Kb, Harzburgite) from the transformation
        //     of the asthenosphere to the basalts-gabbros of the oceanic crust
        //     (p<2.5Kb). Note that the temperature should also be <= LAB_TEMPERATURE_IN_KELVINS
        //     to allow those transformation to take place under an active oceanic ridge.
        static constexpr const double MOHO_PRESSURE_IN_PASCALS= KBARS_2_PASCALS * 2.5; //GPA_2_PASCALS * 0.25;

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
      //  EquationOfState::MulticomponentIncompressible<dim> equation_of_state;
      //
      //  /**
      //   * Object that handles phase transitions.
      //   */
      //  MaterialUtilities::PhaseFunction<dim> phase_function;

    };

  }
}

#endif
