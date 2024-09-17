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

#ifndef PBMFRACTALPARTICLEPHASE_H
#define PBMFRACTALPARTICLEPHASE_H

#include "../particlephase/pbmparticlephase.h"


/// Concept:
/// P --> PBMAggregate

template<typename A>
class PBMFractalParticlePhase : public PBMParticlePhase<A> {

    double D_f; ///< Fractal dimension of the aggregates
    NucleationTheory* nt; ///< Pointer to currently used nucleation theory
    bool init = false; ///< Boolean value which stores whether the model is initialized or not.

public:

    PBMFractalParticlePhase(double _D_f, double _volume = 1e-18)
        : PBMParticlePhase<A>(_volume), D_f(_D_f) {
          // set max and min number of aggregates
          // THIS IS THE DEFAULT. THE USER CAN CHANGE THESE VALUES AFTER THE CLASS INSTANTION
          this->set_max_aggregates(2000);
          this->set_min_aggregates(1990);

          volume = this->get_volume();
        }

    double volume; ///<control volume size

    double get_aggregates_number(){
      return this->aggregates.size();
    }

    void move_aggregates(double dt, double vel, PBMFractalParticlePhase<A> *dest );

private:

    double nucleation(double j);

    void initialize () {
      nt = this-> template findNearest<NucleationTheory>();

      // Mark the object as initialized
      init = true;
    };
};

template<typename A>
double PBMFractalParticlePhase<A>::nucleation(double j) {

        // Initialize the model if not done before
        if (init == false) { initialize(); }

        // create a new aggregate
        std::shared_ptr<Particle> p0(new Particle(j,this->sp,nt));
        std::shared_ptr<A> a0(new A(D_f,p0));

        this->aggregates.push_back(a0);

        return j/this->volume;
}

template<typename A>
void PBMFractalParticlePhase<A>::move_aggregates(double dt, double vel, PBMFractalParticlePhase<A> *dest ) {

  double total_mass = 0.;
  for (auto i : this->aggregates) {
    total_mass += i->get_mass();
  }

  double mass_flux = vel * this->get_volume() * this->get_aggregates_density() * total_mass * dt; //kg
  double flux_mass = 0.; //kg

//  int rem = 0;

  while ( (flux_mass < mass_flux) && (this->aggregates.size() > this->min_aggregates_number) ) {

    int N = this->aggregates.size();

    auto it = this->aggregates.begin();
    double rho = ndm::uniform_double_distr(ndm::rand_gen);
    int idx = int(rho*(N-1));

    std::shared_ptr<A> sel_agg = this->get_aggregate(idx);

    flux_mass += sel_agg->get_mass();

    // Copy the aggregate in the given PBM object if is not a nullptr
    if (dest != nullptr) {
      dest->aggregates.push_back(sel_agg);
      int Nd = dest->aggregates.size();
      dest->volume *= (Nd+1.0)/Nd;
    }

    std::advance(it,idx);

    this->aggregates.erase(it);

    this->volume *= (N-1.0)/N;

//    rem++;

  }

//    std::cout << flux_mass << " of " << mass_flux << " removed." << std::endl;
//    std::cout << rem << " aggregates moved." << std::endl;

  this->aggregates_number_balance();

}

#endif // PBMFRACTALPARTICLEPHASE_H
