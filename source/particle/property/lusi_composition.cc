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


#include <aspect/particle/property/lusi_composition.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      LUSIComposition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                            std::vector<double> &data) const
      {
        Composition<dim>::initialize_one_particle_property(position,data);
      }

      template <int dim>
      void
      LUSIComposition<dim>::update_particle_property(const unsigned int data_position,
                                                 const Vector<double> &solution,
                                                 const std::vector<Tensor<1,dim>> &/*gradients*/,
                                                 typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            const unsigned int solution_component = this->introspection().component_indices.compositional_fields[i];
            //particle->get_properties()[data_position+i] = solution[solution_component];
          }
      }

      template <int dim>
      UpdateTimeFlags
      LUSIComposition<dim>::need_update() const
      {
        return Composition<dim>::need_update(); //update_time_step;
      }

      template <int dim>
      UpdateFlags
      LUSIComposition<dim>::get_needed_update_flags () const
      {
        return Composition<dim>::get_needed_update_flags (); //update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      LUSIComposition<dim>::get_property_information() const
      {

        return Composition<dim>::get_property_information();

        //AssertThrow(this->n_compositional_fields() > 0,
        //            ExcMessage("You have requested the particle property <composition>, "
        //                       "but the number of compositional fields is 0. "
        //                       "Please add compositional fields to your model, or remove "
        //                       "this particle property."));
        //std::vector<std::pair<std::string,unsigned int>> property_information;
        //for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
        //  {
        //    const std::string field_name = this->introspection().name_for_compositional_index(i);
        //    property_information.emplace_back(field_name,1);
        //  }
        //return property_information;
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(LUSIComposition,
                                        "lusi composition",
                                        "Implementation of a plugin in which the particle property is defined by the "
                                        "Laval U. Subductiion Initiation (LUSI, G. Mercier Ph.D. thesis) compositional fields in "
                                        "the model. This can be used to track solid composition transformations and trajectories "
                                        "evolution over time. This is a loan from C. Beaumont's SOPALE model code implementation.")
    }
  }
}
