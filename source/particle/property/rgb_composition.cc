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


#include <aspect/particle/property/rgb_composition.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      void
      RGBComposition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                            std::vector<double> &data) const
      {

	double fields_max_compo= 0.0;
	int dominant_compo_index= -1;

	// --- Subtract 2 from this->n_compositional_fields() to consider only the
	//     relevant material fields (i.e. do not consider the accumulated strain fields
	//     at the end of the fields def. list)
        for (unsigned int i = 0; i < this->n_compositional_fields()-2; i++)
	  {
            double field_compo_value_check= this->get_initial_composition_manager().initial_composition(position,i);
	    
	    if (field_compo_value_check >= fields_max_compo)
	      {
		// --- Update fields_max_compo for the next loop step.
		fields_max_compo= field_compo_value_check;
	        dominant_compo_index= this->introspection().component_indices.compositional_fields[i];	
	      }
	  }
	
          AssertThrow(dominant_compo_index == -1,
                      ExcMessage("Cannot have dominant_compo_index == -1 at this point"));

	  // --- Set the RGB vector of the particle according to the dominant composition:
          particle->get_properties()[data_position]= this->compositions_rgb_codes[dominant_compo_index];
      }

      template <int dim>
      void
      RGBComposition<dim>::update_particle_property(const unsigned int data_position,
                                                    const Vector<double> &solution,
                                                    const std::vector<Tensor<1,dim>> &/*gradients*/,
                                                    typename ParticleHandler<dim>::particle_iterator &particle) const
      {

	unsigned int dominant_composition_index= -1;
	
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            //const unsigned int solution_component = this->introspection().component_indices.compositional_fields[i];
            //particle->get_properties()[data_position+i] = solution[solution_component];
          }
      }

      template <int dim>
      UpdateTimeFlags
      RGBComposition<dim>::need_update() const
      {
        return update_time_step;
      }

      template <int dim>
      UpdateFlags
      RGBComposition<dim>::get_needed_update_flags () const
      {
        return update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      RGBComposition<dim>::get_property_information() const
      {

        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <rgb composition>, "
                               "but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

	// --- Subtract 2 from this->compositions_rgb_codes.size() to get only the
	//     relevant material fields (i.e. do not consider the accumulated strain fields
	//     at the end of the fields def. list)
        AssertThrow(this->n_compositional_fields() == this->compositions_rgb_codes.size()-2,
                    ExcMessage("dimensions mismatch: this->n_compositional_fields() != this->compositions_rgb_codes.size()"));        
	
	const std::vector<std::pair<std::string,unsigned int>> property_information (1,std::make_pair("rgb composition",3));
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
      ASPECT_REGISTER_PARTICLE_PROPERTY(RGBComposition,
                                        "rgb composition",
                                        "Implementation of a plugin in which the particle "
                                        "property is defined by the RGB color define for the dominant compositional field in "
                                        "the model. This can be used to track all the solid compositions "
                                        "evolution over time in the same particles output field")
    }
  }
}
