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

#ifndef _aspect_material_model_equation_of_state_multicomponent_incompressible_si_h
#define _aspect_material_model_equation_of_state_multicomponent_incompressible_si_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/equation_of_state/multicomponent_incompressible.h>
//#include <aspect/material_model/equation_of_state/interface.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace EquationOfState
    {
      using namespace dealii;

      /**
       */
      template <int dim>
      class MulticomponentIncompressibleSI :
          public MaterialModel::EquationOfState::MulticomponentIncompressible<dim>
      {
        public:
          /**
           * A function that computes the output of the equation of state @p out
           * for all compositions and phases, given the inputs in @p in and an
           * index input_index that determines which entry of the vector of inputs is used.
           */
          void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
          //void evaluate(MaterialModel::MaterialModelInputs<dim> &in,
                        const unsigned int input_index,
                        MaterialModel::EquationOfStateOutputs<dim> &out) const;

          /**
           * Return whether the model is compressible or not. Incompressibility
           * does not necessarily imply that the density is constant; rather, it
           * may still depend on temperature or pressure. In the current
           * context, compressibility means whether we should solve the continuity
           * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
           * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
           */
          bool is_compressible () const;

          /**
           * Declare the parameters this class takes through input files.
           * The optional parameter @p default_thermal_expansion determines
           * the default value of the thermal expansivity used in the
           * equation of state.
           */
          static
          void
          declare_parameters (ParameterHandler &prm,
                              const double default_thermal_expansion = 3.5e-5);

          /**
           * Read the parameters this class declares from the parameter file.
           * If @p expected_n_phases_per_composition points to a vector of
           * unsigned integers this is considered the number of phase transitions
           * for each compositional field and will be checked against the parsed
           * parameters.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

        //private:
        //  /**
        //   * Vector of reference densities $\rho_0$ with one entry per composition and phase plus one
        //   * for the background field.
        //   */
        //  std::vector<double> densities;
        //
        //  /**
        //   * The reference temperature $T_0$ used in the computation of the density.
        //   * All components use the same reference temperature.
        //   */
        //  double reference_T;
        //
        //  /**
        //   * Vector of thermal expansivities with one entry per composition and phase plus one
        //   * for the background field.
        //   */
        //  std::vector<double> thermal_expansivities;
        //
        //  /**
        //   * Vector of specific heat capacities with one entry per composition and phase plus one
        //   * for the background field.
        //   */
        //  std::vector<double> specific_heats;
      };
    }
  }
}

#endif
