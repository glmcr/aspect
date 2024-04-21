/*
  Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

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

#include <aspect/particle/property/viscosity.h>
#include <aspect/material_model/visco_plastic.h>
//#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      template <int dim>
      Viscosity<dim>::Viscosity ()
        :
        material_inputs(1,0),
        material_outputs(1,0)
      {}
 
      template <int dim>
      void
      Viscosity<dim>::initialize ()
      {
        material_inputs = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());

        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        AssertThrow((Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model())),
                    ExcMessage("This particle property only makes sense in combination with the viscoelastic or visco_plastic material model."));
      }

      
      template <int dim>
      void
      Viscosity<dim>::initialize_one_particle_property(const Point<dim> &,
                                                       std::vector<double> &data) const
      {
        data.push_back(0.0);
      }


      template <int dim>
      void
      Viscosity<dim>::update_particle_property(const unsigned int data_position,
                                               const Vector<double> & solution,
                                               const std::vector<Tensor<1,dim>> &gradients,
                                               typename ParticleHandler<dim>::particle_iterator &particle) const
      {
	
        material_inputs.position[0] = particle->get_location();

#if DEAL_II_VERSION_GTE(9,4,0)
        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                      &(this->get_dof_handler()));
#else
        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(this->get_triangulation()),
                                                                                      &(this->get_dof_handler()));
#endif
        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
	
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        this->get_material_model().evaluate (material_inputs,material_outputs);

	particle->get_properties()[data_position]= material_outputs.viscosities[0];

      }



      template <int dim>
      UpdateTimeFlags
      Viscosity<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      Viscosity<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      Viscosity<dim>::get_property_information() const
      {
        const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("viscosity",1));
        return property_information;
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(Viscosity,
                                        "viscosity",
                                        "A plugin in which the particle property is defined as the viscosity at the particle current position.")
    }
  }
}
