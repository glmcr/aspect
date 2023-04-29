/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#include <assert.h> 
#include <aspect/particle/utilities_lusi.h>
//#include <aspect/utilities.h>
//#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
//#include <aspect/newton.h>
//#include <aspect/adiabatic_conditions/interface.h>
//#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace ParticleUtilities
    {

      // ---
      ThermodynamicStateMarker::ThermodynamicStateMarker()
      {
        this->moving= false;
      }

      // ---
      ThermodynamicStateMarker::ThermodynamicStateMarker(bool moving)
      {
        this->moving= moving;
      }

      // ---
      bool ThermodynamicStateMarker::isMoving() const
      {
        return this->moving;
      }

      // ---
      bool ThermodynamicStateMarker::isFixed() const
      {
        return not this->moving;
      }

      PTStateMarker::PTStateMarker(): ThermodynamicStateMarker::ThermodynamicStateMarker()
      {
        //ThermodynamicStateMarker::ThermodynamicStateMarker();

         this->temperature= this->pressure= std::numeric_limits<double>::max();
      }

      PTStateMarker::PTStateMarker(double pressure, double temperature): ThermodynamicStateMarker::ThermodynamicStateMarker()
      {
        // ThermodynamicStateMarker::ThermodynamicStateMarker();

        this->pressure= pressure;
        this->temperature= temperature;
      }

      PTStateMarker::PTStateMarker(double pressure, double temperature, bool moving):  ThermodynamicStateMarker::ThermodynamicStateMarker(moving)
      {
        //ThermodynamicStateMarker::ThermodynamicStateMarker(moving);

        //AssertThrow(this->insideValidPressuresRange(pressure)," Invalid pressure -> " + std::to_string(pressure) );
        //AssertThrow(this->insideValidTemperaturesRange(temperature)," Invalid temperature -> " + std::to_string(temperature) );

         this->pressure= pressure;
         this->temperature= temperature;
      }

      // inlined in include file
      //double PTStateMarker::getPressure() const
      //{
      //   return this->pressure;
      //}

      // inlined in include file
      //double PTStateMarker::getTemperature() const
      //{
      //   return this->temperature;
      //}

      bool PTStateMarker::insideValidPressuresRange(double pressure)
      {
        return (pressure < PTStateMarker::MAX_PRESSURE and pressure > PTStateMarker::MIN_PRESSURE);
      }

      bool PTStateMarker::insideValidTemperaturesRange(double temperature)
      {
        return (temperature < PTStateMarker::MAX_TEMPERATURE and temperature > PTStateMarker::MIN_TEMPERATURE);
      }

      bool PTStateMarker::insideValidPressuresRange() const
      {
         return (this->pressure < MAX_PRESSURE
                 and this->pressure > MIN_PRESSURE);
      }

      bool PTStateMarker::insideValidTemperaturesRange() const
      {
         return (this->temperature < MAX_TEMPERATURE
                 and this->temperature > MIN_TEMPERATURE);
      }

      bool PTStateMarker::insideValidPTRanges() const
      {
        return (this->insideValidPressuresRange() and this->insideValidTemperaturesRange());
      }

      ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope()
      {
        this->markersVertices= {};
      }

      // ThermodynamicStateMarkersPolytope::~ThermodynamicStateMarkersPolytope()
      // {
      // 	//for (
      //   //this->markersVertices
      // }      

      // ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(const std::vector<ThermodynamicStateMarker*> tsmVector)
      // //ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(const std::vector<PTStateMarker>& ptsmVector)
      // {
      //   const int nbVertices= tsmVector.size();
      
      // 	assert(nbVertices == 3 || nbVertices == 2);
      
      //   this->markersVertices.clear();
      	 
      // 	if ( nbVertices == 3)
      // 	  {
      // 	     // --- Triangle: Just do a shallow copy of the tsmVector
      //        //     
      // 	    this->markersVertices= tsmVector;
      //   } 
      
      // 	else
      // 	  {
      //     //--- Just two vertices here, build a rectangular polygon with them
      //     //    assuming that the two vertices are forming the diagonal of
      //     //    a rectangle which means that the all the vertices values are
      //     //    reasonably different
      
      //      this->markersVertices.push_back(tsmVector[0]);
      //      this->markersVertices.push_back(tsmVector[1]);

      // 	   // --- Build the other two vertices.
      //      const PTStateMarker ptsm0= static_cast<PTStateMarker>(tsmVector[0]);
	   
      // 	   //this->markersVertices.push_back( PTStateMarker( tsmVector[0].PTStateMarker::getPressure(),
      // 	   //				                   tsmVector[1].PTStateMarker::getTemperature()));
	   
      // 	   //this->markersVertices.push_back( PTStateMarker( ((PTStateMarker)tsmVector[1]).getPressure(),
      // 	   //					          ((PTStateMarker)tsmVector[0]).getTemperature()));
      	    
      // 	  }
      // }
      
      //bool ThermodynamicStateMarkersPolytope::pTAreInside(double pressure, double temperature) const
      //{
      ////AssertThrow(PTStateMarker::insideValidPressuresRange(pressure)," Invalid pressure -> " + std::to_string(pressure) );
      ////AssertThrow(PTStateMarker::insideValidTemperaturesRange(temperature)," Invalid temperature -> " + std::to_string(temperature) );
      //
      // bool isInside= true;
      //
      // // --- TODO: implement code!
      //
      // return isInside;
      //}

    } // --- namespace MaterialUtilities
  }
}

