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

#ifndef AGGREGATE_H
#define AGGREGATE_H

#include "../base/mesoobject.h"
#include "../base/objectcounter.h"
#include "../bond/particlebond.h"
#include "../utilities/ndm_random.h"

#include <list>
#include <memory>
#include <valarray>
#include <map>
#include <cmath>
#include <algorithm>
#include <cfloat>


/// The Aggregate class is a container of connected particles.
/// It provides methods for particles accretion due to surface condensation.
/// Sintering level calculation is provided by the bonds.
/// This type provides also a method for coalescence if sintering level of a bond is
/// beyond a specific limit (e.g. s > 0.95).
///
/// Concepts:
/// P --> Particle

template<typename P>
class Aggregate : public MesoObject, public ObjectCounter<Aggregate<P>> {

protected:

    /// Particles composing the aggregate
    std::list<std::shared_ptr<P>> particles;

    /// Bonds connecting the particles
    std::list<std::shared_ptr<ParticleBond<P>>> bonds;

public:

    /// Constructor
    /// \param p0 pointer to the newly created particle.
    Aggregate(std::shared_ptr<P> p0) { particles.push_back(p0); }

    /// Constructor from a list of particles and bonds
    Aggregate(std::list<std::shared_ptr<P>> _particles,
            std::list<std::shared_ptr<ParticleBond<P>>> _bonds) {
            particles = _particles;
            bonds = _bonds;
    }

    /// Copy constructor
    Aggregate(const Aggregate& a1);

    /// Get aggregate mass as sum of the particles masses [kg]
    double get_mass() final;

    /// Get aggregate volume as sum of the particles volumes [m3]
    double get_volume() final;

    /// Number of particles in the aggregate [#]
    double get_particles_number() const { return particles.size(); }

    /// Get the smallest particle diameter in the aggregate [m]
    double get_particles_smallest_diameter() const;

    /// Get the biggest particle diameter in the aggregate [m]
    double get_particles_biggest_diameter() const;

    /// Get the number of bonds in this aggregate
    int get_bonds_number() const { return bonds.size(); }

    /// Particles mean diameter [m]
    double get_particles_mean_diameter() const;

    /// Aggregate area as sum of the particles surface area [m2]
    double get_surface_area() ;

    /// Get the equivalent spherical surface of the aggregate [m2]
    double get_spherical_surface();

    /// Get the diameter of the equivalent spherical surface [m]
    double get_spherical_diameter() ;

    /// Mean sintering level of the aggregate
    double get_mean_sintering_level() const;

    /// Get particles number accounting for the reduction in size due to sintering
    double get_reduced_particles_number();

    /// Aggregate active area taking into account reduction due to sintering [m2]
    double get_active_surface_area() const;

    /// Get primary particles diameter
    std::valarray<double> get_particles_diameter();

    /// Surface condensation on the particles of this aggregate. (Condensation is
    /// calculated using full particle surface, and is not scaled using the active surface
    /// of the aggregate - TODO)
    /// \param dt timestep size [s]
    /// \param Fs surface condensation rate [#/m2 s]
    double condensation(double dt, double Fs);

    /// Execute a sintering step for each particle in the aggregate.
    /// \param dt timestep [s]
    /// \param T temperature [K]
    void sintering(double dt, double T);

    /// Check if a bond has a sintering level higher than a threshold and provide
    /// coalescence for the particles.
    void coalescence(double max_sint_level);

    ///	Merges two aggregates
    ///	\param std::shared_ptr<Aggregate<P>> a1 aggregate to merge
    ///	\param std::shared_ptr<P> p0 contact particle 1
    ///	\param std::shared_ptr<P> p1 contact particle 2
    void merge(std::shared_ptr<Aggregate<P>> _a1,
               std::shared_ptr<P> _p0, std::shared_ptr<P> _p1);

    ///	Merges two aggregates in a periodic domain
    ///	\param std::shared_ptr<Aggregate<P>> a1 aggregate to merge
    ///	\param std::valarray<double> shift coordinates shift
    ///	\param std::shared_ptr<P> p0 contact particle 1
    ///	\param std::shared_ptr<P> p1 contact particle 2
    void merge_periodic(std::shared_ptr<Aggregate<P>> a1,
                                            std::valarray<double> shift,
                                            std::shared_ptr<P> _p0,
                                            std::shared_ptr<P> _p1);

    void merge(std::shared_ptr<Aggregate<P>> a1);
};

template<typename P>
Aggregate<P>::Aggregate(const Aggregate<P>& a1) : ObjectCounter<Aggregate<P>>(a1) {

  // map p0 and p1
  std::map<std::shared_ptr<P>,std::shared_ptr<P>> hash;

  for(auto &p1: a1.particles) {

    std::shared_ptr<P> p0(new P(*p1));
    particles.push_back(p0);

    hash[p1] = p0;
  }

  for(auto &b1: a1.bonds) {

    std::shared_ptr<P> b1p0 = b1->get_v0();
    std::shared_ptr<P> b1p1 = b1->get_v1();

    double s = b1->get_sintering_level();

    std::shared_ptr<ParticleBond<P>> b0(new ParticleBond<P>(hash[b1p0],hash[b1p1],s));

    bonds.push_back(b0);
  }
}


template<typename P>
double Aggregate<P>::get_mass() {

  double m = 0;

  for(auto &p: particles)
    m += p->get_mass();

  return m;
}


template<typename P>
double Aggregate<P>::get_volume() {

  double V = 0;

  for(auto &p: particles)
    V += p->get_volume();

  return V;
}


template<typename P>
double Aggregate<P>::get_particles_smallest_diameter() const {

  double d;
  double d_min = 1e20;

  for(auto& p: particles) {
    d = p->get_diameter();
    if(d<d_min) d_min = d;
  }

  return d_min;
}

template<typename P>
double Aggregate<P>::get_particles_biggest_diameter() const {

  double d;
  double d_max = DBL_MIN;

  for (auto& p : particles) {
    d = p->get_diameter();
    if (d>d_max) d_max = d;
  }

  return d_max;
}


template<typename P>
double Aggregate<P>::get_particles_mean_diameter() const {

  double d_avg = 0.0;

  for(auto &p: particles)
    d_avg += p->get_diameter();

  return d_avg/particles.size();
}


template<typename P>
double Aggregate<P>::get_surface_area() {

  double A = 0;

  for(auto &p: particles)
    A += p->get_surface();

  return A;
}


template<typename P>
double Aggregate<P>::get_spherical_surface() {

  double V = get_volume();

  return pow(36*M_PI*V*V,1./3.);
}


template<typename P>
double Aggregate<P>::get_spherical_diameter() {

  double S = get_spherical_surface();

  return sqrt(S/M_PI);
}


template<typename P>
double Aggregate<P>::get_mean_sintering_level() const {

  double asl = 0;

  for(auto &b: bonds)
    asl += b->get_sintering_level();

  double n = get_particles_number();

  if(n>1)
    return asl/(n - 1);
  else
    return 0.0;
}


template<typename P>
double Aggregate<P>::get_reduced_particles_number() {

  double S = get_surface_area();
  double V = this->get_volume();

  return S*S*S/(36.0 * M_PI * V*V);
}


template<typename P>
double Aggregate<P>::get_active_surface_area() const {

  double S_sph = get_spherical_surface();
  double s_avg = get_mean_sintering_level();

  double l_avg = pow(1.0/get_particles_number(),1.0/3.0);

  return S_sph/(s_avg*(1-l_avg)+l_avg);
}

template<typename P>
std::valarray<double> Aggregate<P>::get_particles_diameter() {

  std::valarray<double> sizes(this->particles.size());

  int idx = 0;
  for (auto& p :this->particles) {
    sizes[idx] = p->get_diameter();
  }

  return sizes;
}


template<typename P>
double Aggregate<P>::condensation(double dt, double Fs) {

  double n_tot = 0;

  /*
   *   for(auto& p: this->particles) {
   *
   *       // number of molecules condensing on the particle surface
   *       double n = p->get_surface()*Fs*dt;
   *
   *       // n is a floating point number instead of int, assuming that p particle molecules are >> 1
   *       p->add_molecules(n);
   *
   *	n_tot += n;
}
*/
  for (auto p = this->particles.begin(); p != this->particles.end(); ++p) {

    // number of molecules condensing on the particle surface
    double n = (*p)->get_surface()*Fs*dt;

    // n is a floating point number instead of int, assuming that p particle molecules are >> 1
    (*p)->add_molecules(n);

    // if p number of molecules is < 1 erase the particle from the list
    double n_p = (*p)->get_n();
    if (n_p < 0.0) {
      // erase bonds
      for (auto b = this->bonds.begin(); b != this->bonds.end(); ++b) {
        if ((*b)->get_v0()->get_id() == (*p)->get_id() || (*b)->get_v1()->get_id() == (*p)->get_id()) {
          b = this->bonds.erase(b);
        }

        if (b == this->bonds.end())
          break;
      }
      // erase particle
      p = this->particles.erase(p);
    }


    if (p == this->particles.end())
      break;

    n_tot += n;
  }

  return n_tot;
}


template<typename P>
void Aggregate<P>::sintering(double dt, double T) {

  for(auto &b: this->bonds) {

    // temporary parameters (W.J. Menz, M. Kraft, Combustion and Flame 160
    // (2013), 947–958)
    double As = 1.15e13; // pre-exponential factor [s/m4]
    double Es = 27677; // sintering characteristic energy [K]
    double ns = 4; // diameter power

    // take the min diameter (Eggersdorfer)
    double d0 = b->get_v0()->get_diameter();
    double d1 = b->get_v1()->get_diameter();
    double dp = std::min(d0,d1);

    // sintering time formula
    double tau = As*pow(dp,ns)*exp(Es/T);
    //double tau = 1.0e20; // Test Parameter

    // sintering timestep
    b->sintering(dt,tau);
  }

  coalescence(0.99);
}


template<typename P>
void Aggregate<P>::coalescence(double max_sint_level) {

  for(auto it = this->bonds.begin(); it != this->bonds.end(); ++it) {
    if((*it)->get_sintering_level() > max_sint_level) {

      /*std::cout	<< "Particles: " << (*it)->get_v0()->get_id() << "-" << (*it)->get_v1()->get_id()
       *                       << "In Aggregate: " << this->get_id()
       *					<< " Sintering Level: " << (*it)->get_sintering_level() << std::endl
       *					<< " Mass 1:" << (*it)->get_v0()->get_mass()
       *					<< " Mass 2:" << (*it)->get_v1()->get_mass() << std::endl
       *					<< " Diameter 1: " << (*it)->get_v0()->get_diameter()
       *					<< " Diameter 2: " << (*it)->get_v1()->get_diameter()
       *					<< std::endl;
       */
      std::shared_ptr<P> p0 = (*it)->get_v0();
      std::shared_ptr<P> p1 = (*it)->get_v1();

      // if p0 is not the lower index particle, than swap
      if(p0->get_id()>p1->get_id())
        std::swap(p0,p1);

      // add p1 molecules to p0
      p0->add_molecules(p1->get_n());

      // redirect bonds from p1 to p0
      for(auto &bb: this->bonds) {
        if(bb->get_v0() == p1) bb->set_v0(p0);
        if(bb->get_v1() == p1) bb->set_v1(p0);
      }

      // remove the emptied particle p1
      this->particles.remove(p1);

      // trick for deleting the current bond without crash
      if (it == bonds.begin()) {
        this->bonds.remove((*it));
        break;
      } else {
        it = std::prev(it);
        this->bonds.remove((*(std::next(it))));
      }
    }
  }
}


template<typename P>
void Aggregate<P>::merge(std::shared_ptr<Aggregate<P>> a1,
             std::shared_ptr<P> p0, std::shared_ptr<P> p1) {

  // move a1 content inside a0
  this->particles.splice(this->particles.end(),a1->particles);
  bonds.splice(bonds.end(),a1->bonds);

  // create a new bond between p0 and p1
  std::shared_ptr<ParticleBond<P>> b01(new ParticleBond<P>(p0,p1));
  bonds.push_back(b01);

  #ifdef VERBOSE
  std::cout << "NEW AGGREGATE!" << std::endl;
  std::cout << " ID: " << this->get_id() << std::endl;
  std::cout << " Particles: " << std::endl;
  for (auto p = this->particles.begin(); p != this->particles.end(); p++) {
    std::cout << (*p)->get_id() << ",";
  }

  std::cout << std::endl;
  std::cout<<"Bonds:"<< std::endl;
  for (auto b = this->bonds.begin(); b != this->bonds.end(); b++) {
    std::cout << "(" << (*b)->get_v0()->get_id() << "," << (*b)->get_v1()->get_id() << ")" << ",";
  }
  std::cout << std::endl;

  std::cout << "SINTERING LEVEL: " << b01->get_sintering_level() << std::endl;
  std::cout << std::endl;

  std::cout << "RECORDED DISTANCE: " << b01->get_bond_distance() << std::endl;
  std::cout << std::endl;

  #endif

  }
  template<typename P>
  void Aggregate<P>::merge_periodic(std::shared_ptr<Aggregate<P>> a1,
                  std::valarray<double> shift,
          std::shared_ptr<P> p0,
          std::shared_ptr<P> p1) {

    // move a1 content inside a0
    // shift coordinates
    for (auto p = a1->particles.begin(); p != a1->particles.end(); p++) {
      std::valarray<double> x = (*p)->get_x();
      x += shift;
      (*p)->set_x(x);
    }

    // add particles
    this->particles.splice(this->particles.end(), a1->particles);
    bonds.splice(bonds.end(), a1->bonds);

    // create a new bond between p0 and p1
    std::shared_ptr<ParticleBond<P>> b01(new ParticleBond<P>(p0, p1));
    bonds.push_back(b01);

    #ifdef VERBOSE
    std::cout << "NEW AGGREGATE!" << std::endl;
    std::cout << " ID: " << this->get_id() << std::endl;
    std::cout << " Particles: " << std::endl;
    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {
      std::cout << (*p)->get_id() << ",";
    }

    std::cout << std::endl;
    std::cout << "Bonds:" << std::endl;
    for (auto b = this->bonds.begin(); b != this->bonds.end(); b++) {
      std::cout << "(" << (*b)->get_v0()->get_id() << "," << (*b)->get_v1()->get_id() << ")" << ",";
    }
    std::cout << std::endl;

    std::cout << "SINTERING LEVEL: " << b01->get_sintering_level() << std::endl;
    std::cout << std::endl;

    std::cout << "RECORDED DISTANCE: " << b01->get_bond_distance() << std::endl;
    std::cout << std::endl;

    #endif

          }


          template<typename P>
          void Aggregate<P>::merge(std::shared_ptr<Aggregate<P>> a1) {

            double rnd;
            std::shared_ptr<P> p0, p1;

            auto it0 = particles.begin();
            auto it1 = a1->particles.begin();

            rnd = ndm::uniform_double_distr(ndm::rand_gen);
            p0 = *(std::next(it0,int(rnd*(particles.size()-1))));

            do {
              rnd = ndm::uniform_double_distr(ndm::rand_gen);
              p1 = *(std::next(it1,int(rnd*(a1->particles.size()-1))));
            } while (p1==p0);

            merge(a1,p0,p1);
          }

#endif // AGGREGATE_H
