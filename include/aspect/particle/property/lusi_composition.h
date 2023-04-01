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
#include <aspect/simulator_access.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      /**
       * Implementation of a plugin in which the particle property is defined by the
       * Laval U. Subductiion Initiation (LUSI, G. Mercier Ph.D. thesis) compositional fields in
       * the model. This can be used to track solid composition transformations and trajectories
       * evolution over time. This is a loan from C. Beaumont's SOPALE model code implementation.
       *
       * @ingroup ParticleProperties
       */
      template <int dim>
      class LUSIComposition : public Property::Composition<dim> //, public ::aspect::SimulatorAccess<dim>
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
          //static constexpr const char* OLM_ASTH_HYBRID_NID= "olmAsthHybrid";

          // ---
          static constexpr const char* OCEANIC_CRUST_NID= "oceanicCrust";
  
          // ---
          static constexpr const char* OCEANIC_SEDS_NID= "oceanicSeds";
	
          // --- Lithosphere <-> asthenosphere T boundary
          static constexpr const double LAB_TEMPERATURE_IN_KELVINS= 1573.0;

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
	  static constexpr const double SURF_PRESSURE_THRESHOLD_IN_PASCALS= KBARS_2_PASCALS * 0.4;
	
      };
    }
  }
}

#endif
