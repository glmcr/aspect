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


#include <aspect/particle/utilities_lusi.h>
//#include <aspect/initial_composition/interface.h>
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

	// --- Now take care of the ad-hoc material changes
        //     (i.e. rock type transformation depending on
        //     the dynamic and thermodynamic conditions)
        const unsigned int asth_mtl_idx=
          this->introspection().compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID);

        //const unsigned int lith_mtl_idx=
        //  this->introspection().compositional_index_for_name(LITHOSPHERIC_MANTLE_NID);

        const unsigned int ssz_lith_mtl_idx=
          this->introspection().compositional_index_for_name(SSZ_LITHOSPHERIC_MANTLE_NID);
	
        const unsigned int oc_crust_idx=
          this->introspection().compositional_index_for_name(OCEANIC_CRUST_NID);

        const unsigned int ssz_oc_crust_idx=
          this->introspection().compositional_index_for_name(SSZ_OCEANIC_CRUST_NID);	
	
        //const unsigned int olm_asth_hybrid_idx=
        //  this->introspection().compositional_index_for_name(OLM_ASTH_HYBRID_NID);

        const unsigned int oc_seds_idx=
          this->introspection().compositional_index_for_name(OCEANIC_SEDS_NID);

        const unsigned int greenschists_idx=
          this->introspection().compositional_index_for_name(GREENSCHISTS_NID);

        const unsigned int amphibolites_idx=
          this->introspection().compositional_index_for_name(AMPHIBOLITES_NID);	

        const unsigned int granulites_idx=
          this->introspection().compositional_index_for_name(GRANULITES_NID);

        const unsigned int eclogites_idx=
          this->introspection().compositional_index_for_name(ECLOGITES_NID);		

        const unsigned int blueschists_idx=
          this->introspection().compositional_index_for_name(BLUESCHISTS_NID);
	
        const unsigned int coesite_idx=
          this->introspection().compositional_index_for_name(COESITE_NID);

	const double pressure_here= \
	  solution[this->introspection().component_indices.pressure];

	const double temperature_here= \
	  solution[this->introspection().component_indices.temperature];

	//--- pointer shortcut to the particle->get_properties()[data_position]
	//
	//double* const __restrict__ part_compo_props= &particle->get_properties().data()[data_position];
        //double* const part_compo_props= &particle->get_properties().data()[data_position];
	double* const part_compo_props= &particle->get_properties().data()[data_position];

        //if (pressure_here <= MOHO_PRESSURE_IN_PASCALS) {
        // std::cout << "LUSIComposition<dim>::update_particle_property: pressure_here=" << pressure_here << std::endl;
        // std::cout << "LUSIComposition<dim>::update_particle_property: temperature_here=" << temperature_here << std::endl;
        //}

	// --- p,T conditions under which asth. transforms to SSZ crust (from the surface down to the moho)
        if ( asth2SSZCrustPTRect.ptInside(pressure_here,temperature_here))
	  {

            //std::cout << "LUSIComposition<dim>::update_particle_property: oc. seds check surf.: pressure_here="
            //                  << pressure_here << ", temperature_here=" << temperature_here << std::endl;

	   // --- Pour some oc. seds. but only where pressure is < SEDS_POUR_PRESSURE_THRESHOLD_IN_PASCALS
	   if (pressure_here < SEDS_POUR_PRESSURE_THRESHOLD_IN_PASCALS)
	     {
	       part_compo_props[oc_seds_idx] += 0.25; //+= 1.5; //0.75;

	       //--- Keeping oc. seds compo prop between 0.25 and 0.5 here.
               part_compo_props[oc_seds_idx]=
		   std::max(0.25,std::min(0.5,part_compo_props[oc_seds_idx]));

               //std::cout << "LUSIComposition<dim>::update_particle_property: oc. seds check surf. OK: part_compo_props[oc_seds_idx]="
               //              << part_compo_props[oc_seds_idx] << std::endl << std::endl;

               //AssertThrow(false,ExcMessage("LUSIComposition<dim>::update_particle_property: oc. seds check surf.:  Debug stop"));
	     }

	   // --- Transfer particle asth. material (could be 0.0) concentration to
	   //     to the SSZ type of oc. crust.
	   lusiMaterialChange(part_compo_props, asth_mtl_idx, ssz_oc_crust_idx, 0.0, 1.0);

	   //// --- Transfer particle asth. material (could be 0.0) concentration to
	   ////     to the SSZ type of oc. crust.
	   //part_compo_props[ssz_oc_crust_idx] += part_compo_props[asth_mtl_idx];
	   ////--- Keeping compo prop between 0.0 and 1.0
           //part_compo_props[ssz_oc_crust_idx]=
           //   std::max(0.0,std::min(1.0,part_compo_props[ssz_oc_crust_idx]));
	   //// --- Need to set asth. to zero here once its concentration
	   ////     has been transfered to ssz oc. crust.
	   //part_compo_props[asth_mtl_idx]= 0.0;
	  
	  } // --- asth -> ssz oc. crust.

	// --- p,T conditions under which asth. transforms to SSZ oc. lith. mantle (moho to LAB)
	if (asth2SSZOlmPTRect.ptInside(pressure_here,temperature_here))
	  {
            lusiMaterialChange(part_compo_props, asth_mtl_idx, ssz_lith_mtl_idx, 0.0, 1.0);
	  
           // part_compo_props[ssz_lith_mtl_idx] += part_compo_props[asth_mtl_idx];
           // part_compo_props[ssz_lith_mtl_idx]=
           //     std::max(0.0,std::min(1.0,part_compo_props[ssz_lith_mtl_idx]));	   
	   // // --- Need to set asth. to zero here once its concentration
	   // //     has been transfered to ssz oc. lith mantle
	   // part_compo_props[asth_mtl_idx]= 0.0;
	   
	  } // --- // --- asth -> ssz  oc. lith mantle

	// --- p,T conditions under which oc. crust transforms to greenschists facies
	if (greenSchistsPTTri1.ptInside(pressure_here,temperature_here) ||
	    greenSchistsPTTri2.ptInside(pressure_here,temperature_here) )
	  {
           lusiMaterialChange(part_compo_props, oc_crust_idx, greenschists_idx , 0.0, 1.0);
	  
           // part_compo_props[greenschists_idx] += part_compo_props[oc_crust_idx];
           // part_compo_props[greenschists_idx]=
           //     std::max(0.0,std::min(1.0,part_compo_props[greenschists_idx]));	   
	   // // --- Need to set asth. to zero here once its concentration
	   // //     has been transfered to ssz oc. lith mantle
	   // part_compo_props[oc_crust_idx]= 0.0; 
	  }

	// --- p,T conditions under which oc. crust and greenschists transform to amphibolites facies
	if (amphibolitesPTTri1.ptInside(pressure_here,temperature_here) ||
	    amphibolitesPTTri2.ptInside(pressure_here,temperature_here) )
	  {
            lusiMaterialChange(part_compo_props, oc_crust_idx, amphibolites_idx , 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, amphibolites_idx , 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, greenschists and amphibolites transform
	//     to granulite facies
	if (granulitesPTTri1.ptInside(pressure_here,temperature_here) ||
	    granulitesPTTri2.ptInside(pressure_here,temperature_here) )
	  {
            lusiMaterialChange(part_compo_props, oc_crust_idx, granulites_idx , 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, granulites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, amphibolites_idx, granulites_idx, 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, and greenschists transform to blueschists facies
	if (blueschistsPTTri1.ptInside(pressure_here,temperature_here))
	  {
            lusiMaterialChange(part_compo_props, oc_crust_idx, blueschists_idx , 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, blueschists_idx , 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, blueschists, greenschists, amphibolites and
	//     granulites transform to eclogites facies
	if (eclogitesPTTri1.ptInside(pressure_here,temperature_here) ||
	    eclogitesPTTri2.ptInside(pressure_here,temperature_here) )
	  {
            lusiMaterialChange(part_compo_props, oc_crust_idx,     eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, blueschists_idx,  eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, amphibolites_idx, eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, granulites_idx,   eclogites_idx, 0.0, 1.0);
	  }

          // --- Prograde only qtz -> coesite transition.
          if (qtz2CoesPTRect.ptInside(pressure_here,temperature_here))
          {
              lusiMaterialChange(part_compo_props, oc_crust_idx, coesite_idx, 0.0, 1.0);
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
