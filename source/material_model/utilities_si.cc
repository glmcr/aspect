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

#include <aspect/material_model/utilities_si.h>
//#include <aspect/utilities.h>
//#include <deal.II/fe/fe_values.h>
#include <deal.II/base/signaling_nan.h>
//#include <aspect/newton.h>
//#include <aspect/adiabatic_conditions/interface.h>
//#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace MaterialModel
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
    ThermodynamicStateMarker::isMoving() const
    {
      return this->moving;
    }

    // ---
    ThermodynamicStateMarker::isFixed() const
    {
      return not this->moving;
    }

    PTStateMarker::PTStateMarker()
    {
       ThermodynamicStateMarker::ThermodynamicStateMarker();

       this->temperature= this->pressure= std::numeric_limits<double>::max();
    }

    PTStateMarker::PTStateMarker()
    {
       ThermodynamicStateMarker::ThermodynamicStateMarker();

       this->temperature= this->pressure= std::numeric_limits<double>::max();
    }

    PTStateMarker::PTStateMarker(double pressure, double temperature)
    {
       ThermodynamicStateMarker::ThermodynamicStateMarker();

       this->pressure= pressure;
       this->temperature= temperature;
    }

    PTStateMarker::PTStateMarker(double pressure, double temperature, bool moving)
    {
       ThermodynamicStateMarker::ThermodynamicStateMarker(moving);

       AssertThrow(this->insideValidPressuresRange(pressure)," Invalid pressure -> " << pressure );
       AssertThrow(this->insideValidTemperaturesRange(temperature)," Invalid temperature -> " << temperature );

       this->pressure= pressure;
       this->temperature= temperature;
    }

    double PTStateMarker::getPressure() const
    {
       return this->pressure;
    }

    double PTStateMarker::getTemperature() const
    {
       return this->temperature;
    }

    static bool PTStateMarker::insideValidPressuresRange(double pressure)
    {
       return (pressure < MAX_PRESSURE and pressure > MIN_PRESSURE);
    }

    static bool PTStateMarker::insideValidTemparaturesRange(double temperature)
    {
       return (temperature < MAX_TEMPERATURE and temperature > MIN_TEMPERATURE);
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
      return (this->insideValidPressuresRange and this->insideValidTemperaturesRange);
    }

    ThermodynamicStateMarkersPolytope::ThermodynamicStateMarkersPolytope()
    {
      this->markersVertices= {};
    }

    bool ThermodynamicStateMarkersPolytope::pTAreInside(double pressure, double temperature) const
    {
       AssertThrow(PTStateMarker::insideValidPressuresRange(pressure)," Invalid pressure -> " << pressure );
       AssertThrow(PTStateMarker::insideValidTemparaturesRange(temperature)," Invalid temperature -> " << temperature );
    }
  }
}

