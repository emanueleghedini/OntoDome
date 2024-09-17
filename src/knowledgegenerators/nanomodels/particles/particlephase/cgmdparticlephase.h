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

#ifndef CGMDPARTICLEPHASE_H
#define CGMDPARTICLEPHASE_H

#include "particlephase.h"

/// Concept:
/// P --> SpatialAggregate

template<typename A>
class CGMDParticlePhase : public ParticlePhase<A> {

public:

  /// Constructor
  CGMDParticlePhase(double _volume) :ParticlePhase<A>(_volume),
  max_aggregates_number(100),
  min_aggregates_number(90) {}

  /// Simulation Timestep
  //virtual double timestep(double dt, const GasPhase& gp, const NucleationTheory& nt) = 0;

protected:

  /// Evaluates the Langevin motion of the particles
  /// \param	double dt: Simulation's Time Step [s]
  /// \param	GasPhase& gp: Gas Phase data
  /// \param	double T: Temperature [K]
  virtual void motion(double dt, double T) = 0;

  /// Coagulation
  virtual void coagulation() = 0;

  /// Maximum number of aggregates. Further nucleations will lead to the removal of
  /// an aggregate and a volume reduction.
  int max_aggregates_number;

  /// Minumum number of aggregates. Further reductions will lead to
  /// the duplication of a randomly chosen aggregate.
  int min_aggregates_number;

public:

  /// Sets the maximum number of aggregates allowed in the simulation
  ///	\param int _agg_max: maximum number of aggregates in the simulation
  void set_max_agg_number(int _agg_max) { max_aggregates_number = _agg_max; }

  /// Sets the minimum number of aggregates allowed in the simulation
  ///	\param int _agg_min: maximum number of aggregates in the simulation
  void set_min_agg_number(int _agg_min) { min_aggregates_number = _agg_min; }

};

#endif // CGMDPARTICLEPHASE_H
