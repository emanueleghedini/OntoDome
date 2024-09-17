/*
    NanoDome - H2020 European Project NanoDome, GA n.646121
    (www.nanodome.eu, https://github.com/nanodome/nanodome)
    e-mail: Emanuele Ghedini, emanuele.ghedini@unibo.it

    Copyright (C) 2018  Alma Mater Studiorum - Università di Bologna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CNT_H
#define CNT_H

#include <cmath>
#include <limits>
#include <exception>
#include <stdexcept>
#include <iostream>

#include "nucleation.h"

/// Classical Nucleation Theory (CNT) implementation. This class is the implementation of the
/// CNT whose rates are calculated using the properties of the condensing species.
class ClassicalNucleationTheory : public NucleationTheory {

protected:
  double s_mass; ///< Species mass [kg]
  double s_bulk_liq; ///< Species liquid bulk density [kg/m3]
  SingleComponentComposition* species; ///< Pointer to the species for which the CNT instance is defined
  SurfaceTension* sten; ///< Pointer to SurfaceTension object related to species
  SaturationPressure* ssat; ///< Pointer to SurfaceTension object related to species
  SurfaceTensionMaterialRelation* stenm; ///< Pointer to the Material relation used for Surface Tension
  SaturationPressureMaterialRelation* ssatm; ///< Pointer to the Material relation used for Saturation Pressure
  Temperature* T; ///< Pointer to gas temperature
  GasModels* gasmodel; ///< Pointer to currently used GasModel
  bool init = false; ///< Boolean value which stores whether the model is initialized or not.

public:
    ClassicalNucleationTheory() {};

protected:
    /// Model initializer
    void initialize() {
      // Get all the required inputs from relations graph
      gasmodel = findNearest<GasModels>();
      T = gasmodel->findNearest<Temperature>();
      species = findNearest<SingleComponentComposition>();
      sten = species->findNearest<SurfaceTension>();
      ssat = species->findNearest<SaturationPressure>();
      stenm = sten->findNearest<SurfaceTensionMaterialRelation>();
      ssatm = ssat->findNearest<SaturationPressureMaterialRelation>();
      s_mass = *species->findNearest<Mass>()->get_data();
      s_bulk_liq = *species->findNearest<BulkDensityLiquid>()->get_data();

      // Mark the object as initialized
      init = true;

    }

public:
    /// Primary particles formation rate [#/m3 s]
    double nucleation_rate() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Get Temperature value
      double T1 = *T->get_data();

      double rate = 0.;
      //unsigned int fp_control_state = _controlfp(_EM_INEXACT, _MCW_EM);

      // Get the supersaturation ration from GasModel
      double S = gasmodel->get_S(species);

      // check if species is saturated, if not nucleation rate is left to zero
      if(S>1.0) {

          // Update Surface Tension value on species
          stenm->run();
          double s_s_ten = *sten->get_data();

          // Normalized surface tension
          double theta = s_s_ten*get_m_surface()/(K_BOL*T1);
          double ns_sat = gasmodel-> get_n_sat(species);

          double A = (S*ns_sat*ns_sat);
          double B = get_m_volume();
          double C1 = 2.0*s_s_ten;
          double C2 = M_PI*s_mass;
          double D1 = theta - 4.0*pow(theta, 3) / (27.0*pow(log(S), 2));
          double C = sqrt(C1 / C2);
          double D = exp(D1);

          rate = A * B * C * D;

      } else if (S==1) {

          rate = std::numeric_limits<double>::min();
      }

      return rate;
    }

    /// Returns the smallest stable cluster particle number [#]
    double stable_cluster_number() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Get Temperature value
      double T1 = *T->get_data();

      // Default cluser number
      double c_size = 1.0;

      // Get the supersaturation ration from GasPhase
      double S = gasmodel->get_S(species);

      // Update Surface Tension's value on species
      stenm->run();

      // Check if species is saturated; if not, the stable size has no sense and is set to one
      // meaning that the smallest cluster is a single monomer (no-cluster)
      if(S>1) {
          c_size = 2.0 * get_m_surface() * *sten->get_data() / (3*K_BOL*T1*log(S));
          c_size = pow(c_size,3);
      }

      return c_size;
    }

    /// Smallest stable cluster diameter [m]
    double stable_cluster_diameter() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return 2*pow( (3./4.) * get_m_volume() * stable_cluster_number() / M_PI, 1./3.);
    }

    /// Returns the surface condensation rate [#/m2 s]
    double condensation_rate() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Get Temperature value
      double T1 = *T->get_data();

      // Get the supersaturation ration from GasPhase
      double S = gasmodel->get_S(species);

      // Update Saturation Pressure's value on species
      ssatm->run();

      return *ssat->get_data() * (S-1.0) / sqrt(2*M_PI*s_mass*K_BOL*T1);
    }

    /// Returns the molecular Volume [m3]
    double get_m_volume() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return s_mass/(AMU*N_AVO*1000*s_bulk_liq);
    }

    /// Returns the molecular Surface [m2]
    double get_m_surface() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return 4.8360*pow(get_m_volume(),2./3.); }

};

#endif // CNT_H
