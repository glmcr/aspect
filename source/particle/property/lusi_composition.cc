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
                                                 const std::vector<Tensor<1,dim>> &gradients,
                                                 typename ParticleHandler<dim>::particle_iterator &particle) const
      {

	// --- Use the super class method to be sure to update the strain accumulator fields (if any)
	Composition<dim>::update_particle_property(data_position,solution,gradients,particle);

	// --- Now take care of the ad-hoc material changes
        //     (i.e. rock type transformation depending on
        //     the dynamic and thermodynamic conditions)
        const unsigned int asth_mtl_idx=
          this->introspection().compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID);

        const unsigned int lith_mtl_idx=
          this->introspection().compositional_index_for_name(LITHOSPHERIC_MANTLE_NID);

        const unsigned int oc_crust_idx=
          this->introspection().compositional_index_for_name(OCEANIC_CRUST_NID);

        const unsigned int olm_asth_hybrid_idx=
          this->introspection().compositional_index_for_name(OLM_ASTH_HYBRID_NID);

        const unsigned int oc_seds_idx=
          this->introspection().compositional_index_for_name(OCEANIC_SEDS_NID);

	const double pressure_here= \
	  solution[this->introspection().component_indices.pressure];

	const double temperature_here= \
	  solution[this->introspection().component_indices.temperature];

	//--- pointer shortcut to the particle->get_properties()[data_position]
	//
	//double* const __restrict__ part_compo_props= &particle->get_properties().data()[data_position];
        //double* const part_compo_props= &particle->get_properties().data()[data_position];
	double* const part_compo_props= &particle->get_properties().data()[data_position];

        if (temperature_here <= LAB_TEMPERATURE_IN_KELVINS)
          {
            if (pressure_here < SURF_PRESSURE_THRESHOLD_IN_PASCALS)
	      {
		// --- Add oceanic sediments composition to the particles properties when p < SURF_PRESSURE_THRESHOLD_IN_PASCALS
                //     (particle y position must be in a top FE grid cell at such a low pressure) and when the other composition
		//     is oceanic crust (basalts+gabbros) or lithospheric mantle and when its oceanic seds compo is
		//     < 0.75 to ensure that all the top FE cells have a significant proportion of this material to help with the
                //     lubrication at the subduction trench location.

		if ( part_compo_props[oc_crust_idx] > 0.1 ||
		     part_compo_props[lith_mtl_idx] > 0.1 ||
		     part_compo_props[olm_asth_hybrid_idx] > 0.1)
		  {
		    part_compo_props[oc_seds_idx] += 1.5; //0.75;

		    //--- Keeping oc. seds compo prop between 1.0 and 1.5 here.
                    part_compo_props[oc_seds_idx]=
		      std::max(1.0,std::min(1.5,part_compo_props[oc_seds_idx]));
		  }
              }
	    else if (pressure_here <= MOHO_PRESSURE_IN_PASCALS)
	      {
		 //--- asthenosphere and-or hybrid material transforms to basaltic oceanic crust
		 //particle->get_properties()[data_position+oc_crust_idx]=
		 //  particle->get_properties()[data_position+asth_mtl_idx] +
		 //    particle->get_properties()[data_position+olm_asth_hybrid_idx];

		 part_compo_props[oc_crust_idx] +=
	           part_compo_props[asth_mtl_idx] +
		     part_compo_props[olm_asth_hybrid_idx];

		 //--- Keeping compo prop between 0.0 and 1.0
                 part_compo_props[oc_crust_idx]=
		   std::max(0.0,std::min(1.0,part_compo_props[oc_crust_idx]));
		 
		 //--- Set the asthenosphere and hybrid material to 0.0
		 //   since their compositions have been transferred to oceanic crust
		 //particle->get_properties()[data_position+asth_mtl_idx]=
		 //  particle->get_properties()[data_position+olm_asth_hybrid_idx]= 0.0;

		 part_compo_props[asth_mtl_idx]=
		   part_compo_props[olm_asth_hybrid_idx]= 0.0;
		   
	      }
            else if (pressure_here<= OLM_MAX_PRESSURE_IN_PASCALS)
              {
		
		//--- asthenosphere and-or hybrid material transforms to (oceanic or continental)
		//    lithos. mantle
		//particle->get_properties()[data_position+lith_mtl_idx]=
		//  particle->get_properties()[data_position+asth_mtl_idx] +
		//    particle->get_properties()[data_position+olm_asth_hybrid_idx];
		
		 part_compo_props[lith_mtl_idx] +=
	           part_compo_props[asth_mtl_idx] +
		     part_compo_props[olm_asth_hybrid_idx];

		 //--- Keeping compo prop between 0.0 and 1.0
                 part_compo_props[lith_mtl_idx]=
		   std::max(0.0,std::min(1.0,part_compo_props[lith_mtl_idx]));
		
		//--- Set the asthenosphere and hybrid material to 0.0
		//    since their compositions have been transferred to
		//    the lithos. mantle
		//particle->get_properties()[data_position+asth_mtl_idx]=
		//  particle->get_properties()[data_position+olm_asth_hybrid_idx]= 0.0;  
		part_compo_props[asth_mtl_idx]=
		  part_compo_props[olm_asth_hybrid_idx] = 0.0;
		
	      }
	    else 
	      {
		
		//--- asthenosphere transforms to the hybrid material
                //particle->get_properties()[data_position+olm_asth_hybrid_idx]=
		//  particle->get_properties()[data_position+asth_mtl_idx];
                part_compo_props[olm_asth_hybrid_idx] += part_compo_props[asth_mtl_idx];

		 //--- Keeping compo prop between 0.0 and 1.0
                 part_compo_props[olm_asth_hybrid_idx]=
		   std::max(0.0,std::min(1.0,part_compo_props[olm_asth_hybrid_idx]));		
		
		//--- Set the asthenosphere to 0.0 since its composition has been
		//    transferred to the bybrid material.
		//particle->get_properties()[data_position+asth_mtl_idx]= 0.0;  
		part_compo_props[asth_mtl_idx]= 0.0;
		
	      }
	  } //--- if (temperature_here <= LAB_TEMPERATURE_IN_KELVINS)
	
	else
          {
	    //--- Here T > LAB_TEMPERATURE_IN_KELVINS so the hybrid material is
	    //    transforming back to asthenosphere whatever the pressure is.
	    
	    //--- asthenosphere transforms to the hybrid material
            //particle->get_properties()[data_position+asth_mtl_idx]=
	    //  particle->get_properties()[data_position+olm_asth_hybrid_idx];
	    part_compo_props[asth_mtl_idx] +=
	      part_compo_props[olm_asth_hybrid_idx];

	    //--- Keeping compo prop between 0.0 and 1.0
            part_compo_props[asth_mtl_idx]=
	      std::max(0.0,std::min(1.0,part_compo_props[asth_mtl_idx]));		    

	    //--- Set the hybrid material property to 0.0 since its composition has been
	    //    transferred to the asthenospheric material.
	    //particle->get_properties()[data_position+olm_asth_hybrid_idx]= 0.0;
	    part_compo_props[olm_asth_hybrid_idx]= 0.0;  
	  }
	
      } //--- update_particle_property

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

        AssertThrow(this->n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <lusi "
                               "composition>, but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));

        std::vector<std::pair<std::string,unsigned int>> property_information;

        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            std::ostringstream field_name;
            field_name << "lusi " << this->introspection().name_for_compositional_index(i);
            property_information.emplace_back(field_name.str(),1);
          }

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
      ASPECT_REGISTER_PARTICLE_PROPERTY(LUSIComposition,
                                        "lusi composition",
                                        "Implementation of a plugin in which the particle property is defined by the "
                                        "Laval U. Subductiion Initiation (LUSI, G. Mercier Ph.D. thesis) compositional fields in "
                                        "the model. This can be used to track solid composition transformations and trajectories "
                                        "evolution over time. This is a loan from C. Beaumont's SOPALE model code implementation.")
    }
  }
}
