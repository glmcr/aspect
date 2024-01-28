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

#ifndef _aspect_particle_property_lusi_composition_h
#define _aspect_particle_property_lusi_composition_h

//#include <aspect/particle/property/interface.h>
#include <aspect/particle/property/composition.h>
#include <aspect/particle/utilities_lusi.h>
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {

      using namespace Particle::ParticleUtilities;

      /**
       * Implementation of a plugin in which the particle property is defined by the
       * Laval U. Subductiion Initiation (LUSI, G. Mercier Ph.D. thesis) compositional fields in
       * the model. This can be used to track solid composition transformations and trajectories
       * evolution over time. This is a loan from C. Beaumont's SOPALE model code implementation.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class LUSIComposition final : public Property::Composition<dim> //, public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Initialization function. This function is called once at the
           * creation of every particle for every property to initialize its
           * value.
           *
           * @param [in] position The current particle position.
           * @param [in,out] particle_properties The properties of the particle
           * that is initialized within the call of this function. The purpose
           * of this function should be to extend this vector by a number of
           * properties.
           */
          void
          initialize_one_particle_property (const Point<dim> &position,
                                            std::vector<double> &particle_properties) const override;

          /**
           * @copydoc aspect::Particle::Property::Interface::update_particle_property()
           */
          virtual
          void
          update_particle_property (const unsigned int data_position,
                                    const Vector<double> &solution,
                                    const std::vector<Tensor<1,dim>> &gradients,
                                    typename ParticleHandler<dim>::particle_iterator &particle) const override;

          /**
           * This implementation tells the particle manager that
           * we need to update particle properties over time.
           */
          UpdateTimeFlags
          need_update () const override;

          /**
           * Return which data has to be provided to update the property.
           * The pressure and temperature need the values of their variables.
           */
          UpdateFlags
          get_needed_update_flags () const override;

          /**
           * Set up the information about the names and number of components
           * this property requires.
           *
           * @return A vector that contains pairs of the property names and the
           * number of components this property plugin defines.
           */
          std::vector<std::pair<std::string, unsigned int>>
          get_property_information() const override;

	  // ---
          //inline
          static constexpr const char* ASTHENOSPHERIC_MANTLE_NID= "asthenosphere";

          // ---
          //inline
          static constexpr const char* MRB_LITHOSPHERIC_MANTLE_NID= "oceanicLithMantleMRB";

          static constexpr const char* SSZ_LITHOSPHERIC_MANTLE_NID= "oceanicLithMantleSSZ";

          static constexpr const char* SC_LITHOSPHERIC_MANTLE_NID= "SCLM";

          // ---
          static constexpr const char* MRB_OCEANIC_CRUST_NID= "oceanicCrustMRB";

	  static constexpr const char* SSZ_OCEANIC_CRUST_NID= "oceanicCrustSSZ";

          // ---
          static constexpr const char* OCEANIC_SEDS_NID= "oceanicSeds";

          //// --- prograde MTM facies for oc. seds (qtz -> coesite)
          //static constexpr const char* COESITE_NID= "coesite";

	  // --- prograde MTM facies for oc. crust
	  static constexpr const char* GREENSCHISTS_NID= "greenschists";

	  static constexpr const char* AMPHIBOLITES_NID= "amphibolites";

	  static constexpr const char* GRANULITES_NID= "granulites";

          static constexpr const char* ECLOGITES_NID= "eclogites";

	  static constexpr const char* BLUESCHISTS_NID= "blueschists";

	  static constexpr const char* PARTIALLY_MELTED_SSZ_ASTH_NID= "pmeltedSszAsth";
        
          static constexpr const char* PARTIALLY_MELTED_MRB_ASTH_NID= "pmeltedMrbAsth";

          static constexpr const char* ASTH_OLM_HYB_MAT_NID= "asthOLMHybMat";

          static constexpr const char* CONT_UPPER_CRUST_NID= "contUppCrust";
      
          static constexpr const char* CONT_LOWER_CRUST_NID= "contLowrust";

        //static constexpr const char* AMPHIBOLITES_PM_NID= "amphibolitesPM";

          static constexpr const char* ACC_TOTAL_STRAIN_NID= "total_strain";

          static constexpr const char* ACC_NONINIT_PLASTIC_STRAIN_NID= "noninitial_plastic_strain";
        
	  // --- Surf T 273.25
	 static constexpr const double SURF_TEMPERATURE_KELVINS= 273.25;

	  // --- Surf pressure (atmos. pressure at sea level
	  static constexpr const double SURF_ATMOS_PRESSURE_PASCALS= 101500.0;

          // --- Lithosphere <-> asthenosphere T boundary
          // LAB T at 1300C
          static constexpr const double LAB_TEMPERATURE_KELVINS= 1573.0;
          // LAB T at 1250C
          //static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1523.0;
	  // LAB T at 1280C
          //static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1553.0;

	  // --- T at the moho
	  static constexpr const double MOHO_TEMPERATURE_KELVINS= 850.0;

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
          static constexpr const double MOHO_PRESSURE_PASCALS= KBARS_2_PASCALS * 2.5; //GPA_2_PASCALS * 0.25;

          // --- Max. pressure at which oc. lith. mantle can be formed from
          //     the solid residue (i.e. harzburgite) of the partial fusion
          //     of the asth. mantle.
          static constexpr const double OLM_MAX_PRESSURE_PASCALS= KBARS_2_PASCALS * 9.5;

	  // --- Approx lithos. pressure of a 1.5km column of oceanic sediments.
	  static constexpr const double SEDS_POUR_PRESSURE_THRESHOLD_PASCALS= KBARS_2_PASCALS * 0.4;

	  static constexpr const double NO_MTC_ON_DISTANCE_FROM_SIDES= 150e3;

          // --- Temperature under which the cooling asth. transform to an hybrid material
          //     that has the same thermal properties as the OLM but the same WOL rheology
          //     as the asth. (T==1250C). No dependence on pressure for this hybrid material
          static constexpr const double ASTH_OLM_HYBRID_MAT_TEMP_THESHOLD_KELVINS= 1523.0;

          // --- 0.5m/y in m/s
          static constexpr const double ASTH_PARTIAL_MELT_TYPE_VEL_THRESHOLD= 1.5854895991882293e-08;

           // ---  (p,T) tri. zone where asth. in partial SSZ melting state forms
           //      SSZ oc. crust 
	   static const PTStateMarkersTriangle asth2SSZCrustPTTri1;
	   static const PTStateMarkersTriangle asth2SSZCrustPTTri2;
        
           // ---  (p,T) tri. zone where asth. in partial SSZ melting state forms
           //      SSZ oc. lith 
	   static const PTStateMarkersTriangle asth2SSZOlmPTTri1;
	   static const PTStateMarkersTriangle asth2SSZOlmPTTri2;

	   // --- (p,T) tri. zone where asth. is in partial SSZ melting state
	   static const PTStateMarkersTriangle pmSszAsthPTTri1;
	   static const PTStateMarkersTriangle pmSszAsthPTTri2;
	   static const PTStateMarkersTriangle pmSszAsthPTTri3;

           // ---  (p,T) tri. zone where asth. in partial MRB melting state forms
           //      MRB oc. crust         
 	   static const PTStateMarkersTriangle asth2MRBCrustPTTri1;
	   static const PTStateMarkersTriangle asth2MRBCrustPTTri2;
        
           // ---  (p,T) tri. zone where asth. in partial MRB melting state forms
           //      MRB oc. lith
	   static const PTStateMarkersTriangle asth2MRBOlmPTTri1;
	   static const PTStateMarkersTriangle asth2MRBOlmPTTri2;

	   // --- (p,T) tri. zone where asth. is in partial MRB melting state
	   static const PTStateMarkersTriangle pmMrbAsthPTTri1;
	   static const PTStateMarkersTriangle pmMrbAsthPTTri2;
	   static const PTStateMarkersTriangle pmMrbAsthPTTri3;       

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the greenschists facies (prograde only)
	   static const PTStateMarkersTriangle greenSchistsPTTri1;

	   // --- Declare the 2nd p,T triangle where oc. crust material transforms to
	   //     the greenschists facies (prograde only)
	   static const PTStateMarkersTriangle greenSchistsPTTri2;

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the amphibolites facies.
    	   static const PTStateMarkersTriangle amphibolitesPTTri1;

	   // --- Declare the 2nd p,T triangle where oc. crust material transforms to
	   //     the amphibolites facies.
	   static const PTStateMarkersTriangle amphibolitesPTTri2;

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the granulites facies.
	   static const PTStateMarkersTriangle granulitesPTTri1;

	   // --- Declare the 2nd p,T triangle where oc. crust material transforms to
	   //     the granulites facies.
   	   static const PTStateMarkersTriangle granulitesPTTri2;

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
	   static const PTStateMarkersTriangle eclogitesPTTri1;

	   // --- Declare the 2nd p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
	   static const PTStateMarkersTriangle eclogitesPTTri2;
        
	   // --- Declare the 3rd p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
           static const PTStateMarkersTriangle eclogitesPTTri3;

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the blueschists facies.
	   static const PTStateMarkersTriangle blueschistsPTTri1;

	   // --- Declare the 2nd p,T triangle where oc. crust material transforms to
	   //     the blueschists facies.         
           static const PTStateMarkersTriangle blueschistsPTTri2;

	   //// --- Define the 2nd p,T triangle where oc. crust material transforms to
	   ////     the eclogites facies.
	   //const PTStateMarkersTriangle blueschistsPTTri2=
	   //      PTStateMarkersTriangle (PTStateMarker(2.0e9,723.0),PTStateMarker(1.6e9,1273.0),PTStateMarker(2.0e9,1273.0));

	//static const PTStateMarker asth2SSZCrustPTRect00;
	//static const PTStateMarker asth2SSZCrustPTRect11;

	//static const PTStateMarker asth2SSZOlmPTRect00;
	//static const PTStateMarker asth2SSZOlmPTRect11;

        static const PTStateMarker occMtcPT1;
        static const PTStateMarker occMtcPT2;
        static const PTStateMarker occMtcPT3;                                        
        static const PTStateMarker occMtcPT4;
        static const PTStateMarker occMtcPT5;
        static const PTStateMarker occMtcPT6;
        static const PTStateMarker occMtcPT7;
        static const PTStateMarker occMtcPT8;
        static const PTStateMarker occMtcPT9;
        static const PTStateMarker occMtcPT10;
        static const PTStateMarker occMtcPT11;

        static const PTStateMarker sszMtcPT1;
        static const PTStateMarker sszMtcPT2;
        static const PTStateMarker sszMtcPT3;                                        
        static const PTStateMarker sszMtcPT4;
        static const PTStateMarker sszMtcPT5;
        static const PTStateMarker sszMtcPT6;
        static const PTStateMarker sszMtcPT7;
        static const PTStateMarker sszMtcPT8;	
        
      private:
	   static inline void lusiMaterialChange(double* const part_compo_props, int matFromIdx, int matToIdx, double matToMin, double matToMax)
	   {

             // --- NOTE: matToMax MUST be between 0.0 and 1.0
             //     NOTE: matToMin MUST be between 0.0 and 1.0

	     // --- Transfer particle matFrom material (could be 0.0) concentration to
	     //     to the matTo material
	     part_compo_props[matToIdx] += matToMax * part_compo_props[matFromIdx];

	     //--- Keeping matTo compo prop between matToMin and matToMax
             part_compo_props[matToIdx]=
                std::max(matToMin,std::min(matToMax,part_compo_props[matToIdx]));

	     // --- Need to set matFrom to (1.0 - matToMax) * part_compo_props[matFromIdx]
             //     here once its concentration matToMax has been transfered to the destination material.
	     part_compo_props[matFromIdx]= (1.0 - matToMax) * part_compo_props[matFromIdx];  //0.0;
	   }
        
      }; // --- class LUSIComposition

      // --- static PTStateMarker objects definitions for class LUSIComposition for the prograde metam. changes
      //     of the oc. crust.
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT1(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.25e9, 573.0);
       
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT2(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.28e9, 773.0);      

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT3(PTStateMarker::PASCALS_2_MEGA_PASCALS*1.35e9, 773.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT4(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.3e9, 1023.0);     

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT5(PTStateMarker::PASCALS_2_MEGA_PASCALS*1.425e9, 1023.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT6(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.8e9, 573.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT7(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.5e9, 373.0);     

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT8(PTStateMarker::PASCALS_2_MEGA_PASCALS*3.0e9, 473.0);     

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT9(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.3e9,
                                                          LUSIComposition<dim>::LAB_TEMPERATURE_KELVINS);
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT10(PTStateMarker::PASCALS_2_MEGA_PASCALS*1.5e9,
                                                          LUSIComposition<dim>::LAB_TEMPERATURE_KELVINS);
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::occMtcPT11(PTStateMarker::PASCALS_2_MEGA_PASCALS*3.0e9, 773.0);          

      // ---

      // --- Define the 1st p,T triangle object where oc. crust material transforms to
      //     the greenschists facies (prograde only)
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      greenSchistsPTTri1(LUSIComposition<dim>::occMtcPT1, LUSIComposition<dim>::occMtcPT2, LUSIComposition<dim>::occMtcPT3);     

      // --- Define the 2nd p,T triangle object where oc. crust material transforms to
      //     the greenschists facies (prograde only)
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      greenSchistsPTTri2(LUSIComposition<dim>::occMtcPT1, LUSIComposition<dim>::occMtcPT6, LUSIComposition<dim>::occMtcPT3);          


      // --- Define the 1st p,T triangle where oc. crust material transforms to
      //     the amphibolites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      amphibolitesPTTri1(LUSIComposition<dim>::occMtcPT2, LUSIComposition<dim>::occMtcPT4, LUSIComposition<dim>::occMtcPT5);

      // --- Define the 2nd p,T triangle where oc. crust material transforms to
      //     the amphibolites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      amphibolitesPTTri2(LUSIComposition<dim>::occMtcPT2, LUSIComposition<dim>::occMtcPT3, LUSIComposition<dim>::occMtcPT5);

      // --- Define the 1st p,T triangle where oc. crust material transforms to
      //     the granulites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      granulitesPTTri1(LUSIComposition<dim>::occMtcPT4, LUSIComposition<dim>::occMtcPT9,LUSIComposition<dim>::occMtcPT10);      

      // --- Define the 2nd p,T triangle where oc. crust material transforms to
      //     the granulites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      granulitesPTTri2(LUSIComposition<dim>::occMtcPT4, LUSIComposition<dim>::occMtcPT5, LUSIComposition<dim>::occMtcPT10);

      // --- Define the 1st p,T triangle where oc. crust material transforms to
      //     the blueschists facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      blueschistsPTTri1(LUSIComposition<dim>::occMtcPT6, LUSIComposition<dim>::occMtcPT7, LUSIComposition<dim>::occMtcPT8);      

      // --- Define the 2nd p,T triangle where oc. crust material transforms to
      //     the blueschists facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      blueschistsPTTri2(LUSIComposition<dim>::occMtcPT6, LUSIComposition<dim>::occMtcPT3, LUSIComposition<dim>::occMtcPT8);       

      // --- Define the 1st p,T triangle where oc. crust material transforms to
      //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      eclogitesPTTri1(LUSIComposition<dim>::occMtcPT3, LUSIComposition<dim>::occMtcPT8, LUSIComposition<dim>::occMtcPT11);

      // --- Define the 2nd p,T triangle where oc. crust material transforms to
      //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      eclogitesPTTri2(LUSIComposition<dim>::occMtcPT3, LUSIComposition<dim>::occMtcPT5, LUSIComposition<dim>::occMtcPT11);
 
      // --- Define the 3rd p,T triangle where oc. crust material transforms to
      //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      eclogitesPTTri3(LUSIComposition<dim>::occMtcPT5, LUSIComposition<dim>::occMtcPT10, LUSIComposition<dim>::occMtcPT11);
      
      // --- static PTStateMarker objects definitions for class LUSIComposition for:
      //     1. Apparition of partial melts of the upwelling hydrated (SSZ) asthenosphere (pmSSZAsth)
      //     2. The formation of new SSZ oc. lithospheric mantle (sszOLM) with the pmSSZAsth
      //     3. The formation of new SSZ oc. lithospheric mantle (sszOCC) with the pmSSZAsth
      //     
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT1(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::SURF_ATMOS_PRESSURE_PASCALS, 1373.0);
       
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT2(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::SURF_ATMOS_PRESSURE_PASCALS, 2023.0); // --- 1750C
      // 1st value 1500C sszMtcPT2(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::SURF_ATMOS_PRESSURE_PASCALS, 1773.0);      

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT3(PTStateMarker::PASCALS_2_MEGA_PASCALS*5.0e9, 2023.0); // --- 50Kb, 1750C
      // 1st values, 10Kb & 1500C szMtcPT3(PTStateMarker::PASCALS_2_MEGA_PASCALS*1.0e9, 1773.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT4(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::MOHO_PRESSURE_PASCALS, 1403.0); // --- 1130C  
      // 1st value, 1200C sszMtcPT4(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::MOHO_PRESSURE_PASCALS, 1473.0);     

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT5(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::MOHO_PRESSURE_PASCALS, LUSIComposition<dim>::SURF_TEMPERATURE_KELVINS);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT6(PTStateMarker::PASCALS_2_MEGA_PASCALS*LUSIComposition<dim>::SURF_ATMOS_PRESSURE_PASCALS, LUSIComposition<dim>::SURF_TEMPERATURE_KELVINS);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT7(PTStateMarker::PASCALS_2_MEGA_PASCALS*5.0e9, 2023.0);  
      // --- 1st values 5kb, 1300C sszMtcPT7(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.5e9, LUSIComposition<dim>::LAB_TEMPERATURE_KELVINS);     

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::
      sszMtcPT8(PTStateMarker::PASCALS_2_MEGA_PASCALS*5.0e9, LUSIComposition<dim>::SURF_TEMPERATURE_KELVINS); // --- 50Kb
      // 1st value 5kb sszMtcPT8(PTStateMarker::PASCALS_2_MEGA_PASCALS*0.5e9, LUSIComposition<dim>::SURF_TEMPERATURE_KELVINS);

      // --- 1st p,T triangle for pm ssz asth. formation from fast upwelling of hydrated asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      pmSszAsthPTTri1(LUSIComposition<dim>::sszMtcPT1, LUSIComposition<dim>::sszMtcPT2, LUSIComposition<dim>::sszMtcPT4);

      // --- 1st p,T triangle for pm mrb asth. formation from slow upwelling of "dry" asth.
      //     Same as pmSszAsthPTTri1
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::pmMrbAsthPTTri1= pmSszAsthPTTri1;      

      // --- 2nd p,T triangle for pm ssz asth. formation from fast upwelling of hydrated asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      pmSszAsthPTTri2(LUSIComposition<dim>::sszMtcPT2, LUSIComposition<dim>::sszMtcPT4, LUSIComposition<dim>::sszMtcPT7);

      // --- 2nd p,T triangle for pm mrb asth. formation from slow upwelling of dry asth.      
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::pmMrbAsthPTTri2= pmSszAsthPTTri2;      

      // --- 3rd p,T triangle for pm ssz asth. formation from fast upwelling of hydrated asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      pmSszAsthPTTri3(LUSIComposition<dim>::sszMtcPT2, LUSIComposition<dim>::sszMtcPT7, LUSIComposition<dim>::sszMtcPT3);      

      // --- 3rd p,T triangle for pm mrb asth. formation from slow upwelling of dry asth. 
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::pmMrbAsthPTTri3= pmSszAsthPTTri3;
      
      // --- 1st p,T triangle for SSZ OLM formation from pm upwelling ssz asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      asth2SSZOlmPTTri1(LUSIComposition<dim>::sszMtcPT4, LUSIComposition<dim>::sszMtcPT7, LUSIComposition<dim>::sszMtcPT8);

      // --- 1st p,T triangle for MRB OLM formation from pm upwelling mrb asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::asth2MRBOlmPTTri1= asth2SSZOlmPTTri1;     

      // --- 2nd p,T triangle for SSZ OLM formation from upwelling pm ssz asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      asth2SSZOlmPTTri2(LUSIComposition<dim>::sszMtcPT4, LUSIComposition<dim>::sszMtcPT5, LUSIComposition<dim>::sszMtcPT8);

      // --- 2nd p,T triangle for MRB OLM formation from pm upwelling mrb asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::asth2MRBOlmPTTri2= asth2SSZOlmPTTri2;      
      
      // --- 1st p,T triangle for SSZ OCC formation from upwelling pm ssz asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      asth2SSZCrustPTTri1(LUSIComposition<dim>::sszMtcPT1, LUSIComposition<dim>::sszMtcPT5, LUSIComposition<dim>::sszMtcPT4);

      // --- 1st p,T triangle for MRB OCC formation from upwelling pm mrb asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::asth2MRBCrustPTTri1= asth2SSZCrustPTTri1;    

      // --- 2nd p,T triangle for SSZ OCC formation from upwelling pm ssz asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      asth2SSZCrustPTTri2(LUSIComposition<dim>::sszMtcPT1, LUSIComposition<dim>::sszMtcPT5, LUSIComposition<dim>::sszMtcPT6);      

      // --- 2nd p,T triangle for MRB OCC formation from upwelling pm mrb asth.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::asth2MRBCrustPTTri2= asth2SSZCrustPTTri2;
      
    }
  }
}

#endif
