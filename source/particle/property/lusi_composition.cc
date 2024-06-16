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
#include <aspect/boundary_velocity/interface.h>
#include <aspect/boundary_velocity/function.h>
#include <aspect/particle/utilities_lusi.h>
#include <aspect/particle/property/lusi_composition.h>
#include <aspect/initial_composition/interface.h>
//#include <aspect/particle/property/viscoplastic_strain_invariants.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      
     template <int dim>
      LUSIComposition<dim>:: LUSIComposition() :material_inputs(1,0){}

      
      template <int dim>
      void
      LUSIComposition<dim>::initialize ()
      {
        AssertThrow(Plugins::plugin_type_matches<const MaterialModel::ViscoPlastic<dim>>
                    (this->get_material_model()),
                    ExcMessage("This initial condition only makes sense in combination "
                               "with the visco_plastic material model."));

        material_inputs = MaterialModel::MaterialModelInputs<dim>(1,this->n_compositional_fields());
      }
      
      template <int dim> void
      LUSIComposition<dim>::get_strain_data_update(struct strain_data &strain_data_update,
						  const Vector<double> &solution,
						  const std::vector<Tensor<1,dim>> &gradients,
						  typename ParticleHandler<dim>::particle_iterator &particle) const
      {

	// Current timestep
        const double dt = this->get_timestep();

        // Velocity gradients
        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];
        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];
        material_inputs.position[0] = particle->get_location();

        // Calculate strain rate from velocity gradients
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        // Put compositional fields into single variable
        for (unsigned int i = 0; i < this->n_compositional_fields(); i++)
          {
            material_inputs.composition[0][i] = solution[this->introspection().component_indices.compositional_fields[i]];
          }

        // Find out plastic yielding by calling function in material model.
        const MaterialModel::ViscoPlastic<dim> &viscoplastic
          = Plugins::get_plugin_as_type<const MaterialModel::ViscoPlastic<dim>>(this->get_material_model());

        //const bool plastic_yielding = viscoplastic.is_yielding(material_inputs);
	strain_data_update.plastic_yielding= viscoplastic.is_yielding(material_inputs);

	// Calculate strain rate second invariant
        const double edot_ii = std::sqrt(std::max(-second_invariant(deviator(material_inputs.strain_rate[0])), 0.));

        // Calculate strain invariant magnitude over the last time step
        const double strain_update = dt*edot_ii;
	
        if (this->introspection().compositional_name_exists("plastic_strain") && strain_data_update.plastic_yielding == true)
          strain_data_update.plastic_strain = strain_update;

        if (this->introspection().compositional_name_exists("viscous_strain") && strain_data_update.plastic_yielding == false)
	  strain_data_update.viscous_strain = strain_update;

	if (this->introspection().compositional_name_exists("total_strain"))
           strain_data_update.total_strain= strain_update;

	if (this->introspection().compositional_name_exists("noninitial_plastic_strain") && strain_data_update.plastic_yielding == true)
          strain_data_update.noninitial_plastic_strain= strain_update;
      }
      
      // template <int dim>
      // void
      // LUSIComposition<dim>::initialize ()
      // {
      // 	//this->n_components = 0;
      // 	this->ViscoPlasticStrainInvariant<dim>::initialize ();
      // }

      template <int dim>
      void
      LUSIComposition<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                            std::vector<double> &data) const
      {
        Composition<dim>::initialize_one_particle_property(position,data);

        //ViscoPlasticStrainInvariant<dim>::initialize_one_particle_property(position,data);
      }

      template <int dim>
      void
      LUSIComposition<dim>::update_particle_property(const unsigned int data_position,
                                                     const Vector<double> &solution,
                                                     const std::vector<Tensor<1,dim>> &gradients,
                                                     typename ParticleHandler<dim>::particle_iterator &particle) const
      {

        //std::vector<Tensor<1,dim>> dummy; 

        //Composition<dim>::update_particle_property(data_position, solution, gradients, particle);

	//// --- Example of access to another particle property (somehwat a hack but it works)
	//const double simAgeInYears= this->get_time()/year_in_seconds;
	//const unsigned int ageInYearsIdx= this->get_particle_world()
        //  .get_property_manager().get_data_info().get_position_by_field_name("ageInYears");

	// --- Now take care of the ad-hoc material changes
        //     (i.e. rock type transformation depending on
        //     the dynamic and thermodynamic conditions)
        const unsigned int asth_mtl_idx=
	  this->Composition<dim>::introspection().compositional_index_for_name(ASTHENOSPHERIC_MANTLE_NID);

	// --- partially melted ssz asth.
        const unsigned int pm_ssz_asth_mtl_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(PARTIALLY_MELTED_SSZ_ASTH_NID);

	// --- partially melted mrb asth.
        const unsigned int pm_mrb_asth_mtl_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(PARTIALLY_MELTED_MRB_ASTH_NID);

	// ---
        const unsigned int mrb_lith_mtl_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(MRB_LITHOSPHERIC_MANTLE_NID);

        const unsigned int ssz_lith_mtl_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(SSZ_LITHOSPHERIC_MANTLE_NID);

        const unsigned int mrb_oc_crust_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(MRB_OCEANIC_CRUST_NID);

        const unsigned int ssz_oc_crust_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(SSZ_OCEANIC_CRUST_NID);

        const unsigned int oc_seds_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(OCEANIC_SEDS_NID);

        // const unsigned int cont_upp_crust_idx=
        //   this->introspection().compositional_index_for_name(CONT_UPPER_CRUST_NID);

	// const unsigned int cont_low_crust_idx=
        //   this->introspection().compositional_index_for_name(CONT_LOWER_CRUST_NID);       
	
        // const unsigned int sc_lith_mtl_idx=
	//   this->introspection().compositional_index_for_name(SC_LITHOSPHERIC_MANTLE_NID);
	
        const unsigned int greenschists_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(GREENSCHISTS_NID);

        const unsigned int amphibolites_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(AMPHIBOLITES_NID);

        //const unsigned int amphibolitesPM_idx=
        //  this->introspection().compositional_index_for_name(AMPHIBOLITES_PM_NID);

        const unsigned int granulites_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(GRANULITES_NID);

        const unsigned int eclogites_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(ECLOGITES_NID);

        const unsigned int blueschists_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(BLUESCHISTS_NID);

        const unsigned int asth_olm_hyb_mat_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(ASTH_OLM_HYB_MAT_NID);

	// --- NOTE: Assuming here that acc_tot_strain_idx is < acc_ninit_plastic_strain_idx
	//           AND that acc_ninit_plastic_strain_idx = acc_tot_strain_idx + 1
	//           (which means that the acc_ninit_plastic_strain value is located
	//            im memory just after the acc_tot_strain value in the particle data memory space 
        const unsigned int acc_tot_strain_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(ACC_TOTAL_STRAIN_NID);

        const unsigned int acc_ninit_plastic_strain_idx=
          this->Composition<dim>::introspection().compositional_index_for_name(ACC_NONINIT_PLASTIC_STRAIN_NID);
	
	//--- pointer shortcut to the particle->get_properties()[data_position]
	//    which allows to index the values inside it (not clean, but it works)
	double* const part_compo_props= &particle->get_properties().data()[data_position];

	//part_compo_props[ageInYearsIdx]= 0;
	
	//ViscoPlasticStrainInvariant<dim>::update_particle_property(0, solution, gradients, particle);
	//ViscoPlasticStrainInvariant<dim>::update_particle_property(, solution, gradients, particle);

	struct strain_data strain_data_update= {false, 0.0, 0.0, 0.0, 0.0 };

	LUSIComposition<dim>::get_strain_data_update(strain_data_update, solution, gradients, particle);

        part_compo_props[acc_tot_strain_idx] += strain_data_update.total_strain;
	
	if (strain_data_update.plastic_yielding) {
          part_compo_props[acc_ninit_plastic_strain_idx] += strain_data_update.noninitial_plastic_strain;
	}

	const double pressureInPascals_here= \
	  solution[this->Composition<dim>::introspection().component_indices.pressure];

        const double pressureInMPa_here= PTStateMarker::PASCALS_2_MEGA_PASCALS*pressureInPascals_here;

        // --- Kelvins
	const double temperature_here= \
	  solution[this->Composition<dim>::introspection().component_indices.temperature];

	//--- pointer shortcut to the particle->get_properties()[data_position]
	//
	//double* const __restrict__ part_compo_props= &particle->get_properties().data()[data_position];
        //double* const part_compo_props= &particle->get_properties().data()[data_position];
	//double* const part_compo_props= &particle->get_properties().data()[data_position];

        //if (pressureInPascals_here <= MOHO_PRESSURE_IN_PASCALS) {
        // std::cout << "LUSIComposition<dim>::update_particle_property: pressure_here=" << pressure_here << std::endl;
        // std::cout << "LUSIComposition<dim>::update_particle_property: temperature_here=" << temperature_here << std::endl;
        //}

        const types::boundary_id top_num_id=
	  this->Composition<dim>::get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

	//const Tensor<1,dim> bnd_velos=
	//   Function::boundary_velocity(0,particle->get_location());
	//bool extension_stage= false; //bnd_velos[0]
	
        bool in_a_top_cell= false;

#if DEAL_II_VERSION_GTE(9,4,0)
        typename DoFHandler<dim>::active_cell_iterator current_cell=
	    typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),&(this->Composition<dim>::get_dof_handler()));
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
        
	// --- Pour some oc. seds. but only if the marker is in a top cell
	//     and we have MRB or SSZ oc. crust > 0.1 (oceanic setting)
	if (in_a_top_cell && (part_compo_props[mrb_oc_crust_idx] > 0.1 || part_compo_props[ssz_oc_crust_idx] > 0.1 ) ) 
	  {
	    // --- Ensure to always have oc. seds composition at 0.55 in the top (surface) cells
            part_compo_props[oc_seds_idx]= 0.55; //std::min(1.2, std::max(1.2, part_compo_props[oc_seds_idx]));

	    //// --- Need to keep the acc. strains at 0.0 for oc. seds at the surface in the extension stage
	    ////     NOTE: Comment those two next lines for the convergence and slab rollback stage.
	    //part_compo_props[acc_tot_strain_idx]=
	    //  part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
	    
            // --- Limit all other compos between 0.0 and 0.45 
            part_compo_props[mrb_oc_crust_idx]= 
             std::max(0.0,std::min(0.45,part_compo_props[mrb_oc_crust_idx]));

            part_compo_props[ssz_oc_crust_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[ssz_oc_crust_idx]));

            part_compo_props[pm_ssz_asth_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[pm_ssz_asth_mtl_idx]));

	    part_compo_props[pm_mrb_asth_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[pm_mrb_asth_mtl_idx]));

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

            part_compo_props[mrb_lith_mtl_idx]=
             std::max(0.0,std::min(0.45,part_compo_props[mrb_lith_mtl_idx]));

	    part_compo_props[asth_olm_hyb_mat_idx]=
	     std::max(0.0,std::min(0.45,part_compo_props[asth_olm_hyb_mat_idx]));
	}

        // // --- Need to keep the acc. strains at 0.0 for MORB crust in the extension stage.
	// //     NOTE: Comment this code block for the convergence and slab rollback stage.
	// if (part_compo_props[mrb_oc_crust_idx] > 0.5)
	//   {
	//     part_compo_props[acc_tot_strain_idx]=
	//       part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
	//   }
	
	//--- Now check if the marker distance from the sides is far enough
	//    to allow metam. changes because it seems that we have some unwanted
	//    significant pressure oscillations near the sides at distance that are
	//    less than NO_MTC_ON_DISTANCE_FROM_SIDES from them
	const double xPositionMeters= particle->get_location()[0];

        double gridXExtent= NO_MTC_ON_DISTANCE_FROM_SIDES;
	double gridYExtent= NO_MTC_ON_DISTANCE_FROM_SIDES;

        if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->Composition<dim>::get_geometry_model())) {

           const GeometryModel::Box<dim> &box_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->Composition<dim>::get_geometry_model());

           gridXExtent= box_geometry_model.get_extents()[0];
	   gridYExtent= box_geometry_model.get_extents()[1];

        } else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->Composition<dim>::get_geometry_model())) {

           const GeometryModel::TwoMergedBoxes<dim> &two_merged_boxes_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>> (this->Composition<dim>::get_geometry_model());

           gridXExtent= two_merged_boxes_geometry_model.get_extents()[0];
	   gridYExtent= two_merged_boxes_geometry_model.get_extents()[1];

        } else {
            AssertThrow (false,
                         ExcMessage ("geometry model not valid for the lusi composition!"));
        }

        if (gridXExtent <= NO_MTC_ON_DISTANCE_FROM_SIDES) {
           AssertThrow (false,ExcMessage ("Cannot have gridXExtent <= NO_MTC_ON_DISTANCE_FROM_SIDES at this point!"));
        }

	// --- Define the yPositionMeters as being at 7km depth to be sure to
	//     be able to correctly determine if the simulation is at the kinematically imposed extension stage
        const double yPositionMeters= gridYExtent - 7000.0 ; //particle->get_location()[1];
	
	const BoundaryVelocity::Function<dim> & bndFunctionObj=
	  this->get_boundary_velocity_manager().template get_matching_boundary_velocity_model<BoundaryVelocity::Function<dim>>();

	const Point<dim> rigth_side_surf_point(gridXExtent-NO_MTC_ON_DISTANCE_FROM_SIDES,yPositionMeters);
	const Point<dim> left_side_surf_point(NO_MTC_ON_DISTANCE_FROM_SIDES,yPositionMeters);

	const Tensor<1,dim> rigth_bnd_velos= bndFunctionObj.boundary_velocity(0,rigth_side_surf_point);
	const Tensor<1,dim> left_bnd_velos= bndFunctionObj.boundary_velocity(0,left_side_surf_point);

	//const Tensor<1,dim> bnd_velos= bndFunctionObj.boundary_velocity(0,particle->get_location());

        bool in_extension_stage= (rigth_bnd_velos[0] > 0.0 || left_bnd_velos[0] < 0.0) ? true : false;
	
	// bool extension_stage= false;
	// if (xPositionMeters < gridXExtent/2.0) {
	//   // --- Left side of the domain box, x velo should be negative for the extension stage
	//   extension_stage= (bnd_velos[0] < 0.0);	     
	// } else {
	//   // --- Right side of the domain box, x velo should be positive for the extension stage
        //   extension_stage= (bnd_velos[0] > 0.0);
	// }
	
	// --- Need to keep the acc. strains at 0.0 for MORB crust and oc. seds for the extension stage only
	if (part_compo_props[mrb_oc_crust_idx] > 0.5 && in_extension_stage)
	  {
	    part_compo_props[acc_tot_strain_idx]=
	      part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
	  }
	
	if (part_compo_props[oc_seds_idx] > 0.5 && in_extension_stage) {
	  part_compo_props[acc_tot_strain_idx]=
	    part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
	}

        // --- Do the check for the safe distance from the sides:
	if (xPositionMeters <= NO_MTC_ON_DISTANCE_FROM_SIDES ||
            xPositionMeters >= (gridXExtent - NO_MTC_ON_DISTANCE_FROM_SIDES)) {
	    //xPositionMeters >= (box_geometry_model.get_extents()[0] - NO_MTC_ON_DISTANCE_FROM_SIDES)) {

	  // --- Do not do any metam. changes here, just returns.
	  return;
	}

	// // --- Get the vertical velocity at the marker position
	// const double vertical_velo= solution[this->Composition<dim>::introspection().component_indices.velocities[dim-1]];
	
	// // --- Determine if the vertical velo allows the pm asth. of ssz type.
	// const bool pm_asth_ssz_vvelo_ok= (vertical_velo > ASTH_PARTIAL_MELT_SSZ_TYPE_VEL_THRESHOLD) ? true: false;

	//const bool pm_asth_ssz_type= extension_stage ? false: true;
	//const bool pm_asth_ssz_type= ((!extension_stage) && pm_asth_ssz_vvelo_ok) ? true : false;
	const bool pm_asth_ssz_type= !in_extension_stage;

	// // --- Determine which type of pm asth. we have depending on the vertical velo. value 
	// const bool pm_asth_mrb_vvelo_ok= (vertical_velo > ASTH_PARTIAL_MELT_MRB_TYPE_VEL_THRESHOLD) ? true: false;
	
	// //const bool pm_asth_mrb_type= extension_stage ? true : false;
	// const bool pm_asth_mrb_type= ( extension_stage && pm_asth_mrb_vvelo_ok) ? true : false;
	const bool pm_asth_mrb_type= in_extension_stage;
	
	// --- (p,T) and upwelling conditions for which the upwelling hydrated asth. and the hyb. asth. mat.
	//     transforms to partially melted SSZ asthenosphere 
        if ( (pmSszAsthPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
              pmSszAsthPTTri2.ptInside(pressureInMPa_here,temperature_here) ||
              pmSszAsthPTTri3.ptInside(pressureInMPa_here,temperature_here) ||
	      pmSszAsthPTTriMain.ptInside(pressureInMPa_here,temperature_here)) && pm_asth_ssz_type )
	  {
	    lusiMaterialChange(part_compo_props, asth_mtl_idx, pm_ssz_asth_mtl_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, asth_olm_hyb_mat_idx, pm_ssz_asth_mtl_idx, 0.0, 1.0);

	    // --- also transform the partially melted MORB asthenosphere (if any) to partially melted SSZ asthenosphere
	    //     (implies hydratation of the partially melted MORB asthenosphere )
	    lusiMaterialChange(part_compo_props, pm_mrb_asth_mtl_idx,  pm_ssz_asth_mtl_idx, 0.0, 1.0);
	  }

	// --- (p,T) and upwelling conditions for which the upwelling "dry" asth. and the hyb. asth. mat.
	//     transforms to partially melted MORB asthenosphere (no MRB melt in the convergence context)
        if ( (pmMrbAsthPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
              pmMrbAsthPTTri2.ptInside(pressureInMPa_here,temperature_here) ||
              pmMrbAsthPTTri3.ptInside(pressureInMPa_here,temperature_here) ) && pm_asth_mrb_type)
	  {
	    lusiMaterialChange(part_compo_props, asth_mtl_idx, pm_mrb_asth_mtl_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, asth_olm_hyb_mat_idx, pm_mrb_asth_mtl_idx, 0.0, 1.0);
	  }	
	
	// --- (p,T) conditions for which upwelling SSZ asth. partial melts transforms to SSZ crust.
        if ( (asth2SSZCrustPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	      asth2SSZCrustPTTri2.ptInside(pressureInMPa_here,temperature_here)) && !in_extension_stage)
	  {

           const double previousSSZMatContent= part_compo_props[ssz_oc_crust_idx];

	   lusiMaterialChange(part_compo_props, pm_ssz_asth_mtl_idx, ssz_oc_crust_idx, 0.0, 1.0);

           if ( previousSSZMatContent < 0.5 && part_compo_props[ssz_oc_crust_idx] > 0.5 ) {
             // --- reset the accumulated strains to zero for this new ssz mat.
             part_compo_props[acc_tot_strain_idx]= 0.0;
             part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
           }
	  } // --- pm asth -> ssz oc. crust.
	
	// --- (p,T) conditions for which upwelling MORB asth. partial melts transforms to MORB crust.
        if ( (asth2MRBCrustPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	      asth2MRBCrustPTTri2.ptInside(pressureInMPa_here,temperature_here) ) && in_extension_stage)
	  {

           const double previousMRBMatContent= part_compo_props[mrb_oc_crust_idx];

	   lusiMaterialChange(part_compo_props, pm_mrb_asth_mtl_idx, mrb_oc_crust_idx, 0.0, 1.0);

           if ( previousMRBMatContent < 0.5 && part_compo_props[mrb_oc_crust_idx] > 0.5 ) {
             // --- reset the accumulated strains to zero for this new mrb mat.
             part_compo_props[acc_tot_strain_idx]= 0.0;
             part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
           }
	  } // --- pm mrb asth -> mrb oc. crust.

	// --- (p,T) conditions for which upwelling partially melted SSZ asth. transforms to SSZ oc. lith. mantle
	if ( (asth2SSZOlmPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	      asth2SSZOlmPTTri2.ptInside(pressureInMPa_here,temperature_here)) && !in_extension_stage)
	  {

            const double previousSSZMatContent= part_compo_props[ssz_lith_mtl_idx];

	    // --- Transfer particle part. melted ssz asth. material (could be 0.0) concentration to
	    //     to the SSZ type of oc. lith. mantle.	    
	    lusiMaterialChange(part_compo_props, pm_ssz_asth_mtl_idx, ssz_lith_mtl_idx, 0.0, 1.0);

            if ( previousSSZMatContent < 0.5 && part_compo_props[ssz_lith_mtl_idx] > 0.5 ) {
              // --- reset the accumulated strains to zero for this new mrb mat.
              part_compo_props[acc_tot_strain_idx]= 0.0;
              part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
            }
	  } // --- pm mrb asth -> mrb  oc. lith mantle
	
	// --- (p,T) conditions for which upwelling partially melted MRB asth. transforms to MRB oc. lith. mantle
	if ((asth2MRBOlmPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	     asth2MRBOlmPTTri2.ptInside(pressureInMPa_here,temperature_here)) && in_extension_stage)
	  {

            const double previousMRBMatContent= part_compo_props[mrb_lith_mtl_idx];

	    // --- Transfer particle part. melted mrb asth. material (could be 0.0) concentration to
	    //     to the MRB type of oc. lith. mantle.	    
	    lusiMaterialChange(part_compo_props, pm_mrb_asth_mtl_idx, mrb_lith_mtl_idx, 0.0, 1.0);

            if ( previousMRBMatContent < 0.5 && part_compo_props[mrb_lith_mtl_idx] > 0.5 ) {
              // --- reset the accumulated strains to zero for this new ssz mat.
              part_compo_props[acc_tot_strain_idx]= 0.0;
              part_compo_props[acc_ninit_plastic_strain_idx]= 0.0;
            }
	  } // --- pm mrb asth -> rmrb  oc. lith mantle

	// --- p,T conditions under which oc. crust transforms to greenschists facies
	if (greenSchistsPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    greenSchistsPTTri2.ptInside(pressureInMPa_here,temperature_here) )
	  {

	   lusiMaterialChange(part_compo_props, oc_seds_idx, greenschists_idx , 0.0, 1.0);
           lusiMaterialChange(part_compo_props, mrb_oc_crust_idx, greenschists_idx , 0.0, 1.0);
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
            lusiMaterialChange(part_compo_props, mrb_oc_crust_idx, amphibolites_idx, 0.0, 1.0);
            lusiMaterialChange(part_compo_props, ssz_oc_crust_idx, amphibolites_idx, 0.0, 1.0);
	    lusiMaterialChange(part_compo_props, greenschists_idx, amphibolites_idx, 0.0, 1.0);
	  }

	// --- p,T conditions under which oc. crust, greenschists and amphibolites transform
	//     to granulite facies
	if (granulitesPTTri1.ptInside(pressureInMPa_here,temperature_here) ||
	    granulitesPTTri2.ptInside(pressureInMPa_here,temperature_here) )
	  {

	    lusiMaterialChange(part_compo_props, oc_seds_idx, granulites_idx , 0.0, 1.0);
            lusiMaterialChange(part_compo_props, mrb_oc_crust_idx, granulites_idx , 0.0, 1.0);
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
            lusiMaterialChange(part_compo_props, mrb_oc_crust_idx, blueschists_idx , 0.0, 1.0);
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
            lusiMaterialChange(part_compo_props, mrb_oc_crust_idx, eclogites_idx, 0.0, 1.0);
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
           // Need to update both of these to send into material model.
        return update_values | update_gradients;

        //        return Composition<dim>::get_needed_update_flags (); //update_values;
      }

      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      LUSIComposition<dim>::get_property_information() const
      {

        AssertThrow(this->Composition<dim>::n_compositional_fields() > 0,
                    ExcMessage("You have requested the particle property <lusi "
                               "composition>, but the number of compositional fields is 0. "
                               "Please add compositional fields to your model, or remove "
                               "this particle property."));
 
        std::vector<std::pair<std::string,unsigned int>> property_information;

        for (unsigned int i = 0; i < this->Composition<dim>::n_compositional_fields(); i++)
          {
            std::ostringstream field_name;
            field_name << "lusi " << this->Composition<dim>::introspection().name_for_compositional_index(i);
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
                                        "Laval U. Subduction Initiation (LUSI, G. Mercier Ph.D. thesis) compositional fields in "
                                        "the model. This can be used to track solid composition transformations and trajectories "
                                        "evolution over time. This is borrowed from C. Beaumont's SOPALE model code implementation.")
     
    }
  }
}
