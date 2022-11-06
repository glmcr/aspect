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

    //template <int dim>
    //bool
    //ViscoPlasticSI<dim>::
    //is_yielding (const double pressure,
    //             const double temperature,
    //             const std::vector<double> &composition,
    //             const SymmetricTensor<2,dim> &strain_rate) const
    //{
    //  return ViscoPlastic<dim>::is_yielding(pressure,temperature, composition, strain_rate) ; //plastic_yielding;
    //}



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

      //std::cout << std::endl <<"Hello world from ViscoPlasticSI::evaluate() !! debug exit 0" << std::endl;
      //std::exit(0);

    }


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
       ViscoPlastic<dim>::parse_parameters(prm);
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

