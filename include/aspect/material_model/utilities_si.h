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

#ifndef _aspect_material_model_utilities_si_h
#define _aspect_material_model_utilities_si_h

#include <aspect/global.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>

namespace aspect
{
  //template <int dim> class SimulatorAccess;

  namespace MaterialModel
  {
    using namespace dealii;

    //template <int dim> struct MaterialModelOutputs;
    //template <int dim> struct EquationOfStateOutputs;

    /**
     * A namespace in which we define utility functions that
     * might be used in many different places in the material
     * model to prevent code duplication.
     */
    namespace MaterialUtilities
    {

      class ThermodynamicStateMarker
      {

        public:

          static constexpr const double UNDEFINED= std::numeric_limits<double>::max();

          //struct thermodynamicStateMarker initialize();
          ThermodynamicStateMarker();

          ThermodynamicStateMarker(double pressure, double temperature);

        private:

          double pressure;    // --- Pascals
          double temperature; // --- Kelvins

          double specificVolume;
          double specificEntropy;

          double specificEnthalpy;
          double specificGibbsEnergy;
          double specificHelmoltzEnergy;
          double specificInternalEnergy;

      }; // --- class ThermodynamicStateMarker

      // ---
      //template <int nbStateVariables, int nbVertices>
      class ThermodynamicStateMarkersPolytope
      {

       public:

         //unsigned int dimension() const;
         //bool inside(const ThermodynamicStateMarker tsm) const;

         ThermodynamicStateMarkersPolytope(unsigned int nbStateVariables, unsigned int nbVertices);

         bool pTAreInside(double pressure, double temperature) const

       private:

         std::vector<ThermodynamicStateMarker> markersVertices;

      }; // --- class ThermodynamicStateMarkersPolytope
    }
  }
}


#endif
