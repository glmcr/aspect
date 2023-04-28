/*
  Copyright (C) 2011 - 2021 by the authors of the ASPECT code.

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

#ifndef _aspect_particles_utilities_lusi_h
#define _aspect_particles_utilities_lusi_h

#include <aspect/global.h>
#include <deal.II/base/point.h>
//#include <deal.II/base/symmetric_tensor.h>
//#include <deal.II/fe/component_mask.h>
//#include <deal.II/base/signaling_nan.h>
//#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  //template <int dim> class SimulatorAccess;

  namespace Particles
  {

    //using namespace dealii;

    //template <int dim> struct MaterialModelOutputs;
    //template <int dim> struct EquationOfStateOutputs;

    /**
     * A namespace in which we define utility functions that
     * might be used in many different places by the particles
     * related code to prevent code duplication.
     */

    namespace ParticlesUtilities
    {

      // ---
      class ThermodynamicStateMarker
      {

        public:

          //static constexpr const double UNDEFINED= std::numeric_limits<double>::max();

          ThermodynamicStateMarker();

          ThermodynamicStateMarker(bool moving);

          bool isMoving() const;

          bool isFixed() const;

        private:

          bool moving;
      };

      class PTStateMarker: public ThermodynamicStateMarker
      {
        public:

           // --- 220 Kelvins
           static constexpr const double MIN_TEMPERATURE= 220.0;

           // --- CMB temperature (approx.)
           static constexpr const double MAX_TEMPERATURE= 4000.0;

           // --- ATM pressure
           static constexpr const double MIN_PRESSURE= 101.5e3;

           // --- CMB pressure (approx.)
           static constexpr const double MAX_PRESSURE= 140.e9;

           PTStateMarker();
           PTStateMarker(double pressure, double temperature);
           PTStateMarker(double pressure, double temperature, bool moving);

           double getPressure() const;
           double getTemperature() const;

           bool insideValidPTRanges() const;

	   bool insideValidPressuresRange() const;
	   bool insideValidTemperaturesRange() const;

	  static bool insideValidPressuresRange(double pressure);
	  static bool insideValidTemperaturesRange(double temperature);

        private:

          double pressure;    // --- Pascals
          double temperature; // --- Kelvins
      };

      //// ---
      // class EnergyStateMarker: public PTStateMarker
      //{
      //
      //  public:
      //
      //    //struct thermodynamicStateMarker initialize();
      //    EnergyStateMarker();
      //
      //    //EnergyStateMarker(double pressure, double temperature);
      //
      //  private:
      //
      //    double specificVolume;
      //    double specificEntropy;
      //
      //    double specificEnthalpy;
      //    double specificGibbsEnergy;
      //    double specificHelmoltzEnergy;
      //    double specificInternalEnergy;
      //
      //}; // --- class ThermodynamicStateMarker

      // ---
      //template <int nbStateVariables, int nbVertices>
      class ThermodynamicStateMarkersPolytope
      {

       public:

         //unsigned int dimension() const;
         //bool inside(const ThermodynamicStateMarker tsm) const;

         ThermodynamicStateMarkersPolytope();
         //ThermodynamicStateMarkersPolytope(unsigned int nbStateVariables, unsigned int nbVertices);
         ThermodynamicStateMarkersPolytope(unsigned int nbVertices);

         bool pTAreInside(double pressure, double temperature) const;

       private:

         std::vector<ThermodynamicStateMarker> markersVertices;

      }; // --- class ThermodynamicStateMarkersPolytope
    } // --- namespace ParticlesUtilities
  } // --- namespace Particles
} // --- namespace aspect

#endif
