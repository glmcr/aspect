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
#include <deal.II/base/signaling_nan.h>

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

        // std::cout << "PTStateMarker::PTStateMarker constr: pressure=" << pressure << std::endl;
        // std::cout << "PTStateMarker::PTStateMarker constr: temperature=" << temperature  << std::endl;

        this->pressure= pressure;
        this->temperature= temperature;

         //std::cout << "PTStateMarker::PTStateMarker constr: this->pressure=" << this->pressure << std::endl;
         //std::cout << "PTStateMarker::PTStateMarker constr: this->temperature=" << this->temperature  << std::endl;
      }

      PTStateMarker::PTStateMarker(double pressure, double temperature, bool moving):  ThermodynamicStateMarker::ThermodynamicStateMarker(moving)
      {
        //ThermodynamicStateMarker::ThermodynamicStateMarker(moving);

        //AssertThrow(this->insideValidPressuresRange(pressure)," Invalid pressure -> " + std::to_string(pressure) );
        //AssertThrow(this->insideValidTemperaturesRange(temperature)," Invalid temperature -> " + std::to_string(temperature) );

         //std::cout << "PTStateMarker::PTStateMarker constr: pressure=" << pressure << std::endl;
         //std::cout << "PTStateMarker::PTStateMarker constr: pressure=" <<  << std::endl;

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
        //this->markersVertices= {};
        this->nbVertices= 0;
      }

      ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(unsigned int nbVertices)
      {
       this->nbVertices= nbVertices;     
      }
 
      PTStateMarkersTriangle::PTStateMarkersTriangle(): 
         ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope()
      {
	
	this->PTSMSRefs[0]=
	  this->PTSMSRefs[1]=
	    this->PTSMSRefs[2]= NULL;
      }

      PTStateMarkersTriangle::PTStateMarkersTriangle(const PTStateMarker& ptsm0, const PTStateMarker& ptsm1,const PTStateMarker& ptsm2):
               ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(3)
      {
        this->PTSMSRefs[0]= &ptsm0;
	this->PTSMSRefs[1]= &ptsm1;
	this->PTSMSRefs[2]= &ptsm2;
      }

      PTStateMarkersRectangle::PTStateMarkersRectangle():
        ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope()
      {
	this->PTSMSRefs[0]=
	  this->PTSMSRefs[1]= NULL;
      }

      PTStateMarkersRectangle::PTStateMarkersRectangle(const PTStateMarker& ptsm00,const PTStateMarker& ptsm11):
          ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope(2)
      {
        this->PTSMSRefs[0]= &ptsm00; // (p0,T0)
	this->PTSMSRefs[1]= &ptsm11; // (p1,T1)
      }

      // ---
      bool PTStateMarkersRectangle::ptInside(double pressure, double temperature) const
      {

        //std::cout << "PTStateMarkersRectangle::ptInside: pressure=" << pressure << std::endl;
        //std::cout << "PTStateMarkersRectangle::ptInside: temperature=" << temperature << std::endl;
        //std::cout << "PTStateMarkersRectangle::ptInside: this->PTSMSRefs[0]->getPressure()=" << 
          //        this->PTSMSRefs[0]->getPressure() << std::endl;
        //std::cout << "PTStateMarkersRectangle::ptInside: this->PTSMSRefs[0]->getTemperature()=" << 
            //      this->PTSMSRefs[0]->getTemperature() << std::endl;
        //std::cout << "PTStateMarkersRectangle::ptInside: this->PTSMSRefs[1]->getPressure()=" <<
          //        this->PTSMSRefs[1]->getPressure() << std::endl;
        //std::cout << "PTStateMarkersRectangle::ptInside: this->PTSMSRefs[1]->getTemperature()=" << 
          //        this->PTSMSRefs[1]->getTemperature() << std::endl;

        const bool ret= ( (pressure > this->PTSMSRefs[0]->getPressure()) && (pressure < this->PTSMSRefs[1]->getPressure()) &&
                 (temperature > this->PTSMSRefs[0]->getTemperature()) && (temperature < this->PTSMSRefs[1]->getTemperature()) );

        //std::cout << "PTStateMarkersRectangle::ptInside: ret=" <<ret << std::endl << std::endl;
        //AssertThrow(false,ExcMessage("PTStateMarkersRectangle::ptInside: Debug exit"));

        return ret;

        //return ( (pressure > this->PTSMSRefs[0]->getPressure()) && (pressure < this->PTSMSRefs[1]->getPressure()) &&
	//	 (temperature > this->PTSMSRefs[0]->getTemperature()) && (temperature < this->PTSMSRefs[1]->getTemperature()) );
	
	//const bool isInside= (pressure > this->PTSMSRefs[0]->getPressure() && pressure < this->PTSMSRefs[1]->getPressure() &&
	//		      temperature > this->PTSMSRefs[0]->getTemperature() && temperature < this->PTSMSRefs[1]->getTemperature() );
	//return isInside;
      }
      
      // ---
      bool PTStateMarkersTriangle::ptInside(double pressure, double temperature) const
      {
	
        bool isInside= false;

	// --- Local PTStateMarker object to use to see
	//     if (pressure,temperature) combo is inside
	//     this PTStateMarkersTriangle object
	const PTStateMarker ptsmCheck(pressure,temperature);

	// python code
	//denom= crossProd2D(vertex1,vertex2) + crossProd2D(vertex2,vertex3) + crossProd2D(vertex3,vertex1)

	// --- Ordering is important
        const double crossProd01= PTStateMarkerCrossProd(*this->PTSMSRefs[0],*this->PTSMSRefs[1]);
	const double crossProd12= PTStateMarkerCrossProd(*this->PTSMSRefs[1],*this->PTSMSRefs[2]);
	const double crossProd20= PTStateMarkerCrossProd(*this->PTSMSRefs[2],*this->PTSMSRefs[0]);
	
	// --- Ordering is important:
        const double denom= crossProd01 + crossProd12 + crossProd20;
	
	//PTStateMarkerCrossProd(*this->PTSMSRefs[1],*this->PTSMSRefs[2]) +
	//PTStateMarkerCrossProd(*this->PTSMSRefs[2],*this->PTSMSRefs[0]);

        if (std::abs(denom) > PTStateMarker::VERY_SMALL_EPSILON) {

	  //python code
	  //v1MinusV2= ( vertex1[0]-vertex2[0], vertex1[1]-vertex2[1])
	  //v3MinusV1= ( vertex3[0]-vertex1[0], vertex3[1]-vertex1[1])
	  //v2MinusV3= ( vertex2[0]-vertex3[0], vertex2[1]-vertex3[1])

	  // --- Beware! Need to create the PTStateMarker object with the order (p,T) and not (T,p) as in the python code.
	  const PTStateMarker ptsm0MinusPtsm1 (this->PTSMSRefs[0]->getPressure() - this->PTSMSRefs[1]->getPressure(),
					       this->PTSMSRefs[0]->getTemperature() - this->PTSMSRefs[1]->getTemperature());

	  const PTStateMarker ptsm2MinusPtsm0 (this->PTSMSRefs[2]->getPressure() - this->PTSMSRefs[0]->getPressure(),
					       this->PTSMSRefs[2]->getTemperature()  - this->PTSMSRefs[0]->getTemperature());

	  const PTStateMarker ptsm1MinusPtsm2 (this->PTSMSRefs[1]->getPressure() - this->PTSMSRefs[2]->getPressure(),
	 				       this->PTSMSRefs[1]->getTemperature() - this->PTSMSRefs[2]->getTemperature());

	  // python code:
	  //weight12= ( crossProd2D(vertex1,vertex2) + crossProd2D(point2D,v1MinusV2) )/denom
	  //weight31= ( crossProd2D(vertex3,vertex1) + crossProd2D(point2D,v3MinusV1) )/denom
	  //weight23= ( crossProd2D(vertex2,vertex3) + crossProd2D(point2D,v2MinusV3) )/denom      	  

	  const double weight01= ( crossProd01 + PTStateMarkerCrossProd(ptsmCheck, ptsm0MinusPtsm1)) / denom;
	  const double weight20= ( crossProd20 + PTStateMarkerCrossProd(ptsmCheck, ptsm2MinusPtsm0)) / denom;
	  const double weight12= ( crossProd12 + PTStateMarkerCrossProd(ptsmCheck, ptsm1MinusPtsm2)) / denom;

	  // python code:
	  /*if ( weight12 > _minus_epsilon or math.ceil(weight12) == 1) \
	   *and ( weight31 > _minus_epsilon or math.ceil(weight31) == 1 ) \
           *and ( weight23 > _minus_epsilon or math.ceil(weight23) == 1 ) : ret = True
	   */
	  
          isInside= ( ( weight01 > PTStateMarker::NGV_VERY_SMALL_EPSILON || std::ceil(weight01) == 1) && 
		        ( weight20 > PTStateMarker::NGV_VERY_SMALL_EPSILON || std::ceil(weight20) == 1) &&
		          ( weight12 > PTStateMarker::NGV_VERY_SMALL_EPSILON || std::ceil(weight12) == 1) );
	  
	}
	
	return isInside;
      }

    } // --- namespace ParticleUtilities
  } // --- namespace Particle
}

