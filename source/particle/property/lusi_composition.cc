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

#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
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

	// --- partially melted ssz asth.
        const unsigned int pm_ssz_asth_mtl_idx=
          this->introspection().compositional_index_for_name(PARTIALLY_MELTED_SSZ_ASTH_NID);

	// ---
        const unsigned int lith_mtl_idx=
          this->introspection().compositional_index_for_name(LITHOSPHERIC_MANTLE_NID);

        const unsigned int ssz_lith_mtl_idx=
          this->introspection().compositional_index_for_name(SSZ_LITHOSPHERIC_MANTLE_NID);

        const unsigned int oc_crust_idx=
          this->introspection().compositional_index_for_name(OCEANIC_CRUST_NID);

        const unsigned int ssz_oc_crust_idx=
          this->introspection().compositional_index_for_name(SSZ_OCEANIC_CRUST_NID);

        const unsigned int oc_seds_idx=
          this->introspection().compositional_index_for_name(OCEANIC_SEDS_NID);

        const unsigned int greenschists_idx=
          this->introspection().compositional_index_for_name(GREENSCHISTS_NID);

        const unsigned int amphibolites_idx=
          this->introspection().compositional_index_for_name(AMPHIBOLITES_NID);

        //const unsigned int amphibolitesPM_idx=
        //  this->introspection().compositional_index_for_name(AMPHIBOLITES_PM_NID);

        const unsigned int granulites_idx=
          this->introspection().compositional_index_for_name(GRANULITES_NID);

        const unsigned int eclogites_idx=
          this->introspection().compositional_index_for_name(ECLOGITES_NID);

        const unsigned int blueschists_idx=
          this->introspection().compositional_index_for_name(BLUESCHISTS_NID);

        const unsigned int asth_olm_hyb_mat_idx=
          this->introspection().compositional_index_for_name(ASTH_OLM_HYB_MAT_NID);

        //const unsigned int coesite_idx=
	//this->introspection().compositional_index_for_name(COESITE_NID);

	const double pressureInPascals_here= \
	  solution[this->introspection().component_indices.pressure];

        const double pressureInMPa_here= PTStateMarker::PASCALS_2_MEGA_PASCALS*pressureInPascals_here;

        // --- Kelvins
	const double temperature_here= \
	  solution[this->introspection().component_indices.temperature];

	//--- pointer shortcut to the particle->get_properties()[data_position]
	//
	//double* const __restrict__ part_compo_props= &particle->get_properties().data()[data_position];
        //double* const part_compo_props= &particle->get_properties().data()[data_position];
	double* const part_compo_props= &particle->get_properties().data()[data_position];

        //if (pressureInPascals_here <= MOHO_PRESSURE_IN_PASCALS) {
        // std::cout << "LUSIComposition<dim>::update_particle_property: pressure_here=" << pressure_here << std::endl;
        // std::cout << "LUSIComposition<dim>::update_particle_property: temperature_here=" << temperature_here << std::endl;
        //}

        const types::boundary_id top_num_id=
	  this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

        bool in_a_top_cell= false;

#if DEAL_II_VERSION_GTE(9,4,0)
        typename DoFHandler<dim>::active_cell_iterator current_cell=
	    typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),&(this->get_dof_handler()));
#else
        typename DoFHandler<dim>::active_cell_iterator current_cell=
	   typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(this->get_triangulation()),&(this->get_dof_handler()));
#endif

	for (const unsigned int face_no : current_cell->face_indices()) {
	  
           if ( current_cell->face(face_no)->at_boundary() &&
	        current_cell->face(face_no)->boundary_id() == top_num_id) {
	          in_a_top_cell = true;
		  break;
	   }
	}
        
	// --- Pour some oc. seds. but only where pressureInPascals is < SEDS_POUR_PRESSURE_THRESHOLD_IN_PASCALS
	//if (pressureInPascals_here < SEDS_POUR_PRESSURE_THRESHOLD_IN_PASCALS)
	if (in_a_top_cell) 
	  {
            // --- Ensure to always have oc. seds composition at 1.2 in the top (surface) cells
            part_compo_props[oc_seds_idx]= 0.55; //std::min(1.2, std::max(1.2, part_compo_props[oc_seds_idx]));

            // --- Limit all other compos between 0.0 and 0.45 
            part_compo_props[oc_crust_idx]= 
             std::max(0.0,std::min(0.45,part_compo_props[oc_crust_idx]));

            part_compo_props[ssz_oc_crust_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[ssz_oc_crust_idx]));

            part_compo_props[pm_ssz_asth_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[pm_ssz_asth_mtl_idx]));

            part_compo_props[ssz_lith_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[ssz_lith_mtl_idx]));

            part_compo_props[greenschists_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[greenschists_idx]));

            part_compo_props[amphibolites_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[amphibolites_idx]));

            part_compo_props[granulites_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[granulites_idx]));

            part_compo_props[eclogites_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[eclogites_idx]));

            part_compo_props[blueschists_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[blueschists_idx]));

            part_compo_props[asth_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[asth_mtl_idx]));

            part_compo_props[lith_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[lith_mtl_idx]));

            part_compo_props[asth_olm_hyb_mat_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[asth_olm_hyb_mat_idx]));
	}

	//--- Now check if the marker distance from the sides is far enough
	//    to allow metam. changes because it seems that we have some unwanted
	//    significant pressure oscillations near the sides at distance that are
	//    less than NO_MTC_ON_DISTANCE_FROM_SIDES from them
	const double xPositionMeters= particle->get_location()[0];

        double gridXExtent= NO_MTC_ON_DISTANCE_FROM_SIDES;

        if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model())) {

           const GeometryModel::Box<dim> &box_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->get_geometry_model());

           gridXExtent= box_geometry_model.get_extents()[0];

        } else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model())) {

           const GeometryModel::TwoMergedBoxes<dim> &two_merged_boxes_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model());

           gridXExtent= two_merged_boxes_geometry_model.get_extents()[0];

        } else {
            AssertThrow (false,
                         ExcMessage ("geometry model not valid for the lusi composition!"));
        }

        if (gridXExtent <= NO_MTC_ON_DISTANCE_FROM_SIDES) {
           AssertThrow (false,ExcMessage ("Cannot have gridXExtent <= NO_MTC_ON_DISTANCE_FROM_SIDES at this point!"));
        }

        // --- Do the check for the safe distance from the sides:
	if (xPositionMeters <= NO_MTC_ON_DISTANCE_FROM_SIDES ||
            xPositionMeters >= (gridXExtent - NO_MTC_ON_DISTANCE_FROM_SIDES)) {
	    //xPositionMeters >= (box_geometry_model.get_extents()[0] - NO_MTC_ON_DISTANCE_FROM_SIDES)) {

	  // --- Do not do any metam. changes here, just returns.
	  return;
	}

	// --- (p,T) conditions for which the upwelling hydrated asth. (hybrid or not) transforms to partially melted
	//     SSZ asthenosphere 
        if ( pmSszAsthPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
             pmSszAsthPTTri2.ptInside(pressureInMPa_here,temperature_here) ||
             pmSszAsthPTTri3.ptInside(pressureInMPa_here,temperature_here) )

	  {
	    lusiMaterialChange(part_compo_props, asth_mtl_idx, pm_ssz_asth_mtl_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, asth_olm_hyb_mat_idx, pm_ssz_asth_mtl_idx, 0.0, 1.0);

	  }

	// --- (p,T) conditions for which upwelling SSZ asth. partial melts transforms to SSZ crust.
        if ( asth2SSZCrustPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	     asth2SSZCrustPTTri2.ptInside(pressureInMPa_here,temperature_here))
	  {

            //std::cout << "LUSIComposition<dim>::update_particle_property: oc. seds check surf.: pressure_here="
            //                  << pressure_here << ", temperature_here=" << temperature_here << std::endl;

	   // --- Transfer particle part. melted ssz asth. material (could be 0.0) concentration to
	   //     to the SSZ type of oc. crust.
	   //lusiMaterialChange(part_compo_props, asth_mtl_idx, ssz_oc_crust_idx, 0.0, 1.0);
	   lusiMaterialChange(part_compo_props, pm_ssz_asth_mtl_idx, ssz_oc_crust_idx, 0.0, 1.0);

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

	// --- (p,T) conditions for which upwelling partially melted SSZ asth. transforms to SSZ oc. lith. mantle (moho to LAB)
	if (asth2SSZOlmPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    asth2SSZOlmPTTri2.ptInside(pressureInMPa_here,temperature_here))
	  {

	    // --- Transfer particle part. melted ssz asth. material (could be 0.0) concentration to
	    //     to the SSZ type of oc. lith. mantle.	    
	    lusiMaterialChange(part_compo_props, pm_ssz_asth_mtl_idx, ssz_lith_mtl_idx, 0.0, 1.0);
            //lusiMaterialChange(part_compo_props, asth_mtl_idx, ssz_lith_mtl_idx, 0.0, 1.0);
	  
           // part_compo_props[ssz_lith_mtl_idx] += part_compo_props[asth_mtl_idx];
           // part_compo_props[ssz_lith_mtl_idx]=
           //     std::max(0.0,std::min(1.0,part_compo_props[ssz_lith_mtl_idx]));	   
	   // // --- Need to set asth. to zero here once its concentration
	   // //     has been transfered to ssz oc. lith mantle
	   // part_compo_props[asth_mtl_idx]= 0.0;
	   
	  } // --- // --- asth -> ssz  oc. lith mantle

	// --- p,T conditions under which oc. crust transforms to greenschists facies
	if (greenSchistsPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    greenSchistsPTTri2.ptInside(pressureInMPa_here,temperature_here) )
	  {

	   lusiMaterialChange(part_compo_props, oc_seds_idx, greenschists_idx , 0.0, 1.0);
           lusiMaterialChange(part_compo_props, oc_crust_idx, greenschists_idx , 0.0, 1.0);
           lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, greenschists_idx , 0.0, 1.0);
	  
           // part_compo_props[greenschists_idx] += part_compo_props[oc_crust_idx];
           // part_compo_props[greenschists_idx]=
           //     std::max(0.0,std::min(1.0,part_compo_props[greenschists_idx]));	   
	   // // --- Need to set asth. to zero here once its concentration
	   // //     has been transfered to ssz oc. lith mantle
	   // part_compo_props[oc_crust_idx]= 0.0; 
	  }

	// --- p,T conditions under which oc. crust and greenschists transform to amphibolites facies
	if (amphibolitesPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    amphibolitesPTTri2.ptInside(pressureInMPa_here,temperature_here) )
	  {
	    lusiMaterialChange(part_compo_props, oc_seds_idx, amphibolites_idx, 0.0, 1.0);
            lusiMaterialChange(part_compo_props, oc_crust_idx, amphibolites_idx, 0.0, 1.0);
            lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, amphibolites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, amphibolites_idx, 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, greenschists and amphibolites transform
	//     to granulite facies
	if (granulitesPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    granulitesPTTri2.ptInside(pressureInMPa_here,temperature_here) )
	  {

	    lusiMaterialChange(part_compo_props, oc_seds_idx, granulites_idx , 0.0, 1.0);
            lusiMaterialChange(part_compo_props, oc_crust_idx, granulites_idx , 0.0, 1.0);
            lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, granulites_idx , 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, granulites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, amphibolites_idx, granulites_idx, 0.0, 1.0);

	    //// --- Parametrization of the amphibolite facies materials partial melting
	    ////    (%5 partial melt, other materials that amphibolite should be at %0.0
	    ////    presence at those (p,T) conditions but we nevertheless convert them to
	    ////    %95 granulites and %5 partial melt just in case)
	    //lusiMaterialChange(part_compo_props, oc_seds_idx, granulites_idx, 0.0, 0.9);
	    //lusiMaterialChange(part_compo_props, oc_seds_idx, amphibolitesPM_idx, 0.0, 0.1);
            //lusiMaterialChange(part_compo_props, oc_crust_idx, granulites_idx, 0.0, 0.9);
	    //lusiMaterialChange(part_compo_props, oc_crust_idx, amphibolitesPM_idx, 0.0, 0.1);
            // lusiMaterialChange(part_compo_props, greenschists_idx, granulites_idx, 0.0, 0.9);
	    //lusiMaterialChange(part_compo_props, greenschists_idx, amphibolitesPM_idx, 0.0, 0.1);

	    //lusiMaterialChange(part_compo_props, amphibolites_idx, granulites_idx, 0.0, 0.9);
	    //lusiMaterialChange(part_compo_props, amphibolites_idx, amphibolitesPM_idx, 0.0, 0.1);
	  }

	// --- p,T conditions under which oc. crust, and greenschists transform to blueschists facies
	if (blueschistsPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
            blueschistsPTTri2.ptInside(pressureInMPa_here,temperature_here))
	  {
	    lusiMaterialChange(part_compo_props, oc_seds_idx, blueschists_idx , 0.0, 1.0);
            lusiMaterialChange(part_compo_props, oc_crust_idx, blueschists_idx , 0.0, 1.0);
            lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, blueschists_idx , 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, blueschists_idx , 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, blueschists, greenschists, amphibolites and
	//     granulites transform to eclogites facies
	if (eclogitesPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    eclogitesPTTri2.ptInside(pressureInMPa_here,temperature_here) ||
            eclogitesPTTri3.ptInside(pressureInMPa_here,temperature_here) )
	  {
	    lusiMaterialChange(part_compo_props, oc_seds_idx,      eclogites_idx, 0.0, 1.0);
            lusiMaterialChange(part_compo_props, oc_crust_idx,     eclogites_idx, 0.0, 1.0);
            lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, blueschists_idx,  eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, amphibolites_idx, eclogites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, granulites_idx,   eclogites_idx, 0.0, 1.0);
	  }

        // --- cooled asth -> asth OLM hybrid.
        if (temperature_here < ASTH_OLM_HYBRID_MAT_TEMP_THESHOLD_KELVINS)
	  {
            lusiMaterialChange(part_compo_props, asth_mtl_idx, asth_olm_hyb_mat_idx, 0.0, 1.0);
          }
	else // --- heated asth OLM hybrid (if any) -> asth
	  {
            lusiMaterialChange(part_compo_props, asth_olm_hyb_mat_idx, asth_mtl_idx, 0.0, 1.0);
	  }

        //// --- Prograde only oc. seds. (i.e. mainly qtz) -> coesite transition.
        //if (qtz2CoesPTRect.ptInside(pressure_here,temperature_here))
        //{
        //   lusiMaterialChange(part_compo_props, oc_seds_idx, coesite_idx, 0.0, 1.0);
        //}

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
