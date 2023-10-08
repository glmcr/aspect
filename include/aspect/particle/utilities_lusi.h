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
//#include <aspect/simulator_access.h>
//#include <deal.II/base/symmetric_tensor.h>
//#include <deal.II/fe/component_mask.h>
//#include <deal.II/base/signaling_nan.h>
//#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  //template <int dim> class SimulatorAccess;

  namespace Particle
  {

     //using namespace dealii;
     //using namespace dealii::Particles;

    //template <int dim> struct MaterialModelOutputs;
    //template <int dim> struct EquationOfStateOutputs;

    /**
     * A namespace in which we define utility functions that
     * might be used in many different places by the particles
     * related code to prevent code duplication.
     */

    namespace ParticleUtilities
    {

      //static inline double crossProd2D(

      // ---
      class ThermodynamicStateMarker
      {

        public:

          //static constexpr const double UNDEFINED= std::numeric_limits<double>::max();

          ThermodynamicStateMarker();

          ThermodynamicStateMarker(bool moving);

          bool isMoving() const;

          bool isFixed() const;

	//virtual double getPressure() const = 0;
	//virtual double getTemperature() const = 0;

        private:

          bool moving;
      };

      class PTStateMarker: public ThermodynamicStateMarker
      {
        public:

	   static constexpr const double VERY_SMALL_EPSILON= 1e-75;

	   static constexpr const double NGV_VERY_SMALL_EPSILON= -VERY_SMALL_EPSILON;

           // --- 220 Kelvins
           static constexpr const double MIN_TEMPERATURE_KELVINS= 220.0;

           // --- CMB temperature (approx.)
           static constexpr const double MAX_TEMPERATURE_KELVINS= 4000.0;

           // --- ATM pressure
           static constexpr const double MIN_PRESSURE_PASCALS= 101.5e3;

           // --- CMB pressure (approx.)
           static constexpr const double MAX_PRESSURE_PASCALS= 140.e9;

           static constexpr const double MEGA_PASCALS_2_PASCALS= 1e6;

           static constexpr const double PASCALS_2_MEGA_PASCALS= 1.0/MEGA_PASCALS_2_PASCALS;

           static constexpr const double GIGA_PASCALS_2_MEGA_PASCALS= 1e3;

           static constexpr const double MEGA_PASCALS_2_GIGA_PASCALS= 1.0/GIGA_PASCALS_2_MEGA_PASCALS;

           PTStateMarker();
           PTStateMarker(double pressureInMPa, double temperatureInK);
           PTStateMarker(double pressureInMPa, double temperatureInK, bool moving);

	   inline virtual double getPressureInMPa() const final { return pressureInMPa; }
	   inline virtual double getTemperatureInK() const final { return temperatureInK; }

           inline virtual double getPressureInPa() const final { return MEGA_PASCALS_2_PASCALS*pressureInMPa; }
           inline virtual double getTemperatureInC() const final { return temperatureInK - 273.25; }

           inline virtual double getPressureInGPa() const final { return MEGA_PASCALS_2_GIGA_PASCALS*pressureInMPa; }

           bool insideValidPTRanges() const;

	   bool insideValidPressuresRange() const;
	   bool insideValidTemperaturesRange() const;

	   static bool insideValidPressuresRange(double pressureInPascals);
	   static bool insideValidTemperaturesRange(double temperatureInK);

        private:

          double pressureInMPa; // --- Mega Pascals
          double temperatureInK; // --- Kelvins
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

         //static constexpr const double PASCALS_2_MEGA_PASCALS= 1.0/10e6;

         //unsigned int dimension() const;
         //bool inside(const ThermodynamicStateMarker tsm) const;

         ThermodynamicStateMarkersPolytope();
	//~ThermodynamicStateMarkersPolytope();
	
         //ThermodynamicStateMarkersPolytope(unsigned int nbStateVariables, unsigned int nbVertices);
         ThermodynamicStateMarkersPolytope(unsigned int nbVertices);

	  //ThermodynamicStateMarkersPolytope(const std::vector<ThermodynamicStateMarker*> tsmVector);
	  //ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(const std::vector<PTStateMarker>&);
	
          //bool pTAreInside(double pressure, double temperature) const;
          //virtual bool ptsmInside(const ThermodynamicStateMarker&) const= 0;
	  virtual bool ptInside(double pressureInMPa, double temperatureInK) const= 0;

         protected:

           unsigned int nbVertices;

	//private:
        // std::vector<ThermodynamicStateMarker> markersVertices;
         

      }; // --- class ThermodynamicStateMarkersPolytope

      // class ThermodynamicStateMarkersTriangle: public  ThermodynamicStateMarkersPolytope {};
      // class ThermodynamicStateMarkersRectangle: public  ThermodynamicStateMarkersPolytope {};

      class PTStateMarkersTriangle: public ThermodynamicStateMarkersPolytope
      {
        public:
	  PTStateMarkersTriangle();
	  PTStateMarkersTriangle(const PTStateMarker& ptsm0, const PTStateMarker& ptsm1, const PTStateMarker& ptsm2);
	
	  virtual bool ptInside(double pressureInMPa, double temperatureInK) const final override;
	
        private:
	  const PTStateMarker* PTSMSRefs [3]; //= { null, null, null };	  
      };

      class PTStateMarkersRectangle: public ThermodynamicStateMarkersPolytope
      {

	//   Regular p,T, rectangle p1 > p0, T1 > T0
	//    (p1,T0) --------- (p1,T1)
	//       |                 | 
	//       |                 |
	//    (p0,T0) --------- (p0,T1)

        public:
	  PTStateMarkersRectangle();
	  PTStateMarkersRectangle(const PTStateMarker& ptsm0, const PTStateMarker& ptsm1);

	  virtual bool ptInside(double pressureInMPa, double temperatureInK) const final override;
	
        private:
	  const PTStateMarker* PTSMSRefs [2];  	
      };
      
      static inline double PTStateMarkerCrossProd(const PTStateMarker& PTSM1, const PTStateMarker& PTSM2) {
	return PTSM1.getTemperatureInK() * PTSM2.getPressureInMPa() - PTSM2.getTemperatureInK() * PTSM1.getPressureInMPa();
      }
	
    } // --- namespace ParticleUtilities
  } // --- namespace Particle
} // --- namespace aspect

#endif
