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
          static constexpr const char* LITHOSPHERIC_MANTLE_NID= "oceanicLithMantle";

          static constexpr const char* SSZ_LITHOSPHERIC_MANTLE_NID= "oceanicLithMantleSSZ";

          //// --- Oceanic lithos.  <-> asthenoshperic hybrid material
          ////     (asth. flow law but OLM thermal cond. to parametrize
          ////     the conversion of convecting asth. mantle to lithospheric
          ////     mantle by cooling only i.e. not by partial fusion)
          ////     It does not exists at the beginning and is produced
          ////     when the asthenospheric material cools because of
          ////     the downward migration of the 1573K (1300C) isotherm
          ////     and the adiabatic heating is positive (implying
          ////     very slow downward material advection which would
          ////     be slower than the isotherm downward migration itself)
          //static constexpr const char* OLM_ASTH_HYBRID_NID= "olmAsthHybrid"

          // ---
          static constexpr const char* OCEANIC_CRUST_NID= "oceanicCrust";

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

	  // --- Surf T 273.25
	  static constexpr const double SURF_TEMPERATURE= 273.25;

	  // --- Surf pressure (atmos. pressure at sea level
	  static constexpr const double SURF_ATMOS_PRESSURE= 101500.0;

          // --- Lithosphere <-> asthenosphere T boundary
          // LAB T at 1300C
          static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1573.0;
          // LAB T at 1250C
          //static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1523.0;
	  // LAB T at 1280C
          //static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1553.0;

	  // --- T at the moho
	  static constexpr const double MOHO_TEMPERATURE_IN_KELVINS= 850.0;

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
          static constexpr const double MOHO_PRESSURE_IN_PASCALS= KBARS_2_PASCALS * 2.5; //GPA_2_PASCALS * 0.25;

          // --- Max. pressure at which oc. lith. mantle can be formed from
          //     the solid residue (i.e. harzburgite) of the partial fusion
          //     of the asth. mantle.
          static constexpr const double OLM_MAX_PRESSURE_IN_PASCALS= KBARS_2_PASCALS * 9.5;

	  // --- Approx lithos. pressure of a 1.5km column of oceanic sediments.
	  static constexpr const double SEDS_POUR_PRESSURE_THRESHOLD_IN_PASCALS= KBARS_2_PASCALS * 0.4;

          ////--- oc. seds. transforms (mainly) to coesite
          ////    (mainly dependant on pressure and not really on T)
          //static constexpr const double QTZ_TO_COESITE_TEMPERATURE_THRESHOLD_IN_KELVINS= SURF_TEMPERATURE;

          //// ---
          //static constexpr const double QTZ_TO_COESITE_PRESSURE_THRESHOLD_IN_PASCALS= GPA_2_PASCALS*2.0;

	//private:

	   // --- Declare the PTStateMarkersRectangle object that defines the p,T conditions where
	   //     asth. transform to oc. crust.
           //static const PTStateMarkersRectangle asth2SSZCrustPTRect;

	   // --- Declare the PTStateMarkersRectangle object that defines the p,T conditions where
	   //     asth. transform to oc. lith. mantle.
	   //static const PTStateMarkersRectangle asth2SSZOlmPTRect;

	   static const PTStateMarkersTriangle asth2SSZCrustPTTri1;
	   static const PTStateMarkersTriangle asth2SSZCrustPTTri2;

	   static const PTStateMarkersTriangle asth2SSZOlmPTTri1;
	   static const PTStateMarkersTriangle asth2SSZOlmPTTri2;

	   // --- (p,T) tri. zone where ssz asth. in partial melting state
	   static const PTStateMarkersTriangle pmSszAsthPTTri;

	   //// ---
           //static const PTStateMarkersRectangle qtz2CoesPTRect;

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

	   // --- Declare the 1st p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
	static const PTStateMarkersTriangle blueschistsPTTri1;

	   //// --- Define the 2nd p,T triangle where oc. crust material transforms to
	   ////     the eclogites facies.
	   //const PTStateMarkersTriangle blueschistsPTTri2=
	   //      PTStateMarkersTriangle (PTStateMarker(2.0e9,723.0),PTStateMarker(1.6e9,1273.0),PTStateMarker(2.0e9,1273.0));

	//static const PTStateMarker asth2SSZCrustPTRect00;
	//static const PTStateMarker asth2SSZCrustPTRect11;

	//static const PTStateMarker asth2SSZOlmPTRect00;
	//static const PTStateMarker asth2SSZOlmPTRect11;

       static const PTStateMarker pmSszAsthPTTri11;
       static const PTStateMarker pmSszAsthPTTri12;
       static const PTStateMarker pmSszAsthPTTri13;

       static const PTStateMarker asth2SSZCrustPTTri11;
       static const PTStateMarker asth2SSZCrustPTTri12;
       static const PTStateMarker asth2SSZCrustPTTri13;

       static const PTStateMarker asth2SSZCrustPTTri21;
       static const PTStateMarker asth2SSZCrustPTTri22;
       static const PTStateMarker asth2SSZCrustPTTri23;

       static const PTStateMarker asth2SSZOlmPTTri11;
       static const PTStateMarker asth2SSZOlmPTTri12;
       static const PTStateMarker asth2SSZOlmPTTri13;

       static const PTStateMarker asth2SSZOlmPTTri21;
       static const PTStateMarker asth2SSZOlmPTTri22;
       static const PTStateMarker asth2SSZOlmPTTri23;

       //static const PTStateMarker qtz2CoesPTRect00;
       //static const PTStateMarker qtz2CoesPTRect11;

       static const PTStateMarker greenSchistsPTTri10;
       static const PTStateMarker greenSchistsPTTri11;
       static const PTStateMarker greenSchistsPTTri12;

       static const PTStateMarker greenSchistsPTTri20;
       static const PTStateMarker greenSchistsPTTri21;
       static const PTStateMarker greenSchistsPTTri22;

       static const PTStateMarker amphibolitesPTTri10;
       static const PTStateMarker amphibolitesPTTri11;
       static const PTStateMarker amphibolitesPTTri12;

       static const PTStateMarker amphibolitesPTTri20;
       static const PTStateMarker amphibolitesPTTri21;
       static const PTStateMarker amphibolitesPTTri22;

       static const PTStateMarker granulitesPTTri10;
       static const PTStateMarker granulitesPTTri11;
       static const PTStateMarker granulitesPTTri12;

       static const PTStateMarker granulitesPTTri20;
       static const PTStateMarker granulitesPTTri21;
       static const PTStateMarker granulitesPTTri22;

       static const PTStateMarker eclogitesPTTri10;
       static const PTStateMarker eclogitesPTTri11;
       static const PTStateMarker eclogitesPTTri12;

       static const PTStateMarker eclogitesPTTri20;
       static const PTStateMarker eclogitesPTTri21;
       static const PTStateMarker eclogitesPTTri22;

       static const PTStateMarker blueschistsPTTri10;
       static const PTStateMarker blueschistsPTTri11;
       static const PTStateMarker blueschistsPTTri12;

      private:
	   static inline void lusiMaterialChange(double* const part_compo_props, int matFromIdx, int matToIdx, double matToMin, double matToMax)
	   {
	     // --- Transfer particle matFrom material (could be 0.0) concentration to
	     //     to the matTo material
	     part_compo_props[matToIdx] += part_compo_props[matFromIdx];

	     //--- Keeping matTo compo prop between matToMin and matToMax
             part_compo_props[matToIdx]=
                std::max(matToMin,std::min(matToMax,part_compo_props[matToIdx]));

	     // --- Need to set matFrom to zero here once its concentration
	     //     has been transfered to ssz oc. crust.
	     part_compo_props[matFromIdx]= 0.0;
	   }
      };

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::pmSszAsthPTTri11(LUSIComposition<dim>::SURF_ATMOS_PRESSURE,1373.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::pmSszAsthPTTri12(LUSIComposition<dim>::SURF_ATMOS_PRESSURE,
								 LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::pmSszAsthPTTri13(1.5e9,LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);
      //const PTStateMarker LUSIComposition<dim>::pmSszAsthPTTri13(0.95e9,LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);

      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      pmSszAsthPTTri(LUSIComposition<dim>::pmSszAsthPTTri11,
	    	     LUSIComposition<dim>::pmSszAsthPTTri12,
		     LUSIComposition<dim>::pmSszAsthPTTri13);

      // --- Tri1 (p,T) points for SSZ crust creation from partial melting of Asth.
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri11(LUSIComposition<dim>::SURF_ATMOS_PRESSURE,
								     LUSIComposition<dim>::SURF_TEMPERATURE);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri12(LUSIComposition<dim>::SURF_ATMOS_PRESSURE,1373.0);
								     //LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri13(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,1403.0);
								     // LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);

      // --- Define the (p,T) 1st triangle where the SSZ oc. crust is formed with
      //     partial melting of the upwelling asthenosphere.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
         asth2SSZCrustPTTri1(LUSIComposition<dim>::asth2SSZCrustPTTri11,
	    		     LUSIComposition<dim>::asth2SSZCrustPTTri12,
			     LUSIComposition<dim>::asth2SSZCrustPTTri13);

     // --- Tri2 (p,T) points for SSZ crust creation with partial melting of Asth.
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri21(LUSIComposition<dim>::SURF_ATMOS_PRESSURE,
								     LUSIComposition<dim>::SURF_TEMPERATURE);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri22(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,
								     LUSIComposition<dim>::SURF_TEMPERATURE);
								     //1423.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZCrustPTTri23(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,1403.0);
								     //LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);

      // --- Define the (p,T) 2nd triangle where the SSZ oc. crust is formed with
      //     partial melting of the upwelling asthenosphere.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
         asth2SSZCrustPTTri2(LUSIComposition<dim>::asth2SSZCrustPTTri21,
	    		     LUSIComposition<dim>::asth2SSZCrustPTTri22,
			     LUSIComposition<dim>::asth2SSZCrustPTTri23);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri11(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,
								   LUSIComposition<dim>::SURF_TEMPERATURE); //1423.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri12(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,1403.0);
                                                                   //LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri13(0.95e9, 1473.0);
								   //LUSIComposition<dim>::SURF_TEMPERATURE);

      // --- Define the 2 (p,T) triangles where the SSZ oc. lith mantle is formed with
      //     the partial melted ssz upwelling asthenosphere.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
         asth2SSZOlmPTTri1(LUSIComposition<dim>::asth2SSZOlmPTTri11,
	    		   LUSIComposition<dim>::asth2SSZOlmPTTri12,
			   LUSIComposition<dim>::asth2SSZOlmPTTri13);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri21(LUSIComposition<dim>::MOHO_PRESSURE_IN_PASCALS,
								   LUSIComposition<dim>::SURF_TEMPERATURE); //1423.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri22(0.95e9,
								   LUSIComposition<dim>::SURF_TEMPERATURE);
                                                                   //LUSIComposition<dim>::LAB_TEMPERATURE_IN_KELVINS);
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::asth2SSZOlmPTTri23(0.95e9, 1473.0);

      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
         asth2SSZOlmPTTri2(LUSIComposition<dim>::asth2SSZOlmPTTri21,
	    		   LUSIComposition<dim>::asth2SSZOlmPTTri22,
			   LUSIComposition<dim>::asth2SSZOlmPTTri23);

      //// ---- Qtz (oc. seds) to coesite section, (p,T) rectangle
      //template <int dim>
      //const PTStateMarker LUSIComposition<dim>::
      //   qtz2CoesPTRect00(LUSIComposition<dim>::QTZ_TO_COESITE_PRESSURE_THRESHOLD_IN_PASCALS,
      //                    LUSIComposition<dim>::QTZ_TO_COESITE_TEMPERATURE_THRESHOLD_IN_KELVINS);

    //template <int dim>
    //const PTStateMarker LUSIComposition<dim>::
    //   qtz2CoesPTRect11(LUSIComposition<dim>::GPA_2_PASCALS*8.0,2000.0);

    //template <int dim>
    //const PTStateMarkersRectangle LUSIComposition<dim>::
    //    qtz2CoesPTRect(LUSIComposition<dim>::qtz2CoesPTRect00,LUSIComposition<dim>::qtz2CoesPTRect11);

      // --- Greenschists section
      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri10(0.15e9,573.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri11(0.2e9,773.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri12(1.2e9,773.0);
      //const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri12(0.95e9,773.0);

      // --- Define the 1st p,T triangle object where oc. crust material transforms to
      //     the greenschists facies (prograde only)
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
          greenSchistsPTTri1(LUSIComposition<dim>::greenSchistsPTTri10,
			     LUSIComposition<dim>::greenSchistsPTTri11,
			     LUSIComposition<dim>::greenSchistsPTTri12);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri20(0.15e9,573.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri21(0.6e9, 600.0);

      template <int dim>
      const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri22(1.2e9,773.0);
      //const PTStateMarker LUSIComposition<dim>::greenSchistsPTTri22(0.95e9,773.0);
      
       // --- Define the 2nd p,T triangle object where oc. crust material transforms to
       //     the greenschists facies (prograde only)
       template <int dim>
       const PTStateMarkersTriangle LUSIComposition<dim>::
          greenSchistsPTTri2(LUSIComposition<dim>::greenSchistsPTTri20,
	    		     LUSIComposition<dim>::greenSchistsPTTri21,
			     LUSIComposition<dim>::greenSchistsPTTri22);      

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri10(0.2e9, 773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri11(0.25e9,973.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri12(1.25e9,973.0);
      
       // --- Define the 1st p,T triangle where oc. crust material transforms to
       //     the amphibolites facies.      
       template <int dim>
       const PTStateMarkersTriangle LUSIComposition<dim>::
          amphibolitesPTTri1(LUSIComposition<dim>::amphibolitesPTTri10,
	                     LUSIComposition<dim>::amphibolitesPTTri11,
	   		     LUSIComposition<dim>::amphibolitesPTTri12);
       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri20(0.2e9,773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri21(1.2e9,773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::amphibolitesPTTri22(1.25e9,973.0);
      
       // --- Define the 2nd p,T triangle where oc. crust material transforms to
       //     the amphibolites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      amphibolitesPTTri2(LUSIComposition<dim>::amphibolitesPTTri20,
			 LUSIComposition<dim>::amphibolitesPTTri21,
			 LUSIComposition<dim>::amphibolitesPTTri22);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri10(0.25e9,973.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri11(1.25e9,973.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri12(1.4e9,1573.0);
      
      
      // --- Define the 1st p,T triangle where oc. crust material transforms to
      //     the granulites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
            granulitesPTTri1(LUSIComposition<dim>::granulitesPTTri10,
			     LUSIComposition<dim>::granulitesPTTri11,
			     LUSIComposition<dim>::granulitesPTTri12);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri20(0.25e9,973.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri21(0.2e9,1573.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::granulitesPTTri22(1.4e9,1573.0);
      
      // --- Define the 2nd p,T triangle where oc. crust material transforms to
      //     the granulites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
            granulitesPTTri2(LUSIComposition<dim>::granulitesPTTri20,
			     LUSIComposition<dim>::granulitesPTTri21,
			     LUSIComposition<dim>::granulitesPTTri22);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri10(1.2e9,773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri11(3.0e9,673.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri12(3.0e9,1573.0);
      
	   // --- Define the 1st p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      eclogitesPTTri1(LUSIComposition<dim>::eclogitesPTTri10,
		      LUSIComposition<dim>::eclogitesPTTri11,
		      LUSIComposition<dim>::eclogitesPTTri12);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri20(1.2e9,773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri21(3.0e9,1573.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::eclogitesPTTri22(1.4e9,1573.0);
      
	   // --- Define the 2nd p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      eclogitesPTTri2(LUSIComposition<dim>::eclogitesPTTri20,
		      LUSIComposition<dim>::eclogitesPTTri21,
		      LUSIComposition<dim>::eclogitesPTTri22);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::blueschistsPTTri10(1.2e9,773.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::blueschistsPTTri11(0.8e9,473.0);

       template <int dim>
       const PTStateMarker LUSIComposition<dim>::blueschistsPTTri12(3.0e9,673.0);

      
	   // --- Define the 1st p,T triangle where oc. crust material transforms to
	   //     the eclogites facies.
      template <int dim>
      const PTStateMarkersTriangle LUSIComposition<dim>::
      blueschistsPTTri1(LUSIComposition<dim>::blueschistsPTTri10,
			LUSIComposition<dim>::blueschistsPTTri11,
			LUSIComposition<dim>::blueschistsPTTri12);
    }
  }
}

#endif
