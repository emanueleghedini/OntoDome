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

#ifndef PARTICLE_H
#define PARTICLE_H

#include "../base/mesoobject.h"
#include "../base/objectcounter.h"

#include <memory>


/// Particle implements the data structure needed to describe a primary particle
/// and the basic functions to change its composition and to describe it.
/// A particle is a collection of molecules of different species. The
/// particle shape is assumed to be spherical.
class Particle : public MesoObject, public ObjectCounter<Particle> {

public:

    double n;  ///< number of monomers
    SingleComponentComposition* s; ///< species composing the particle
    NucleationTheory* nt;

    /// Constructor.
    /// \param _n number of molecules in the particle [#]
    /// \param _nt pointer to nucleation theory
    /// \param _s pointer to species type
    Particle(double _n, SingleComponentComposition* _s, NucleationTheory* _nt) : n(_n), s(_s), nt(_nt) {

      nt = s->findNearest<NucleationTheory>();

    }


    /// Copy constructor
    Particle( Particle& p1);

    /// Get particle molecules number [#].
    double get_n()  { return n; }

    /// Get particle species.
    SingleComponentComposition* get_species()  { return s; }

    /// Add molecules.
    /// \param dn amount of molecules to add [#]
    void add_molecules(double dn) { n += dn; }

    /// Get particle mass from the molecules composition [kg]
    double get_mass() ;

    /// Get particle volume using the molecules volume [m3]
    double get_volume() ;

    /// Get particle diameter [m]
    double get_diameter() ;

    /// Get particle surface (spherical assumption) [m2]
    double get_surface() ;
};

Particle::Particle(Particle& p1) : ObjectCounter<Particle>(p1), s(p1.s), nt(p1.nt) {

    n = p1.n;
}


double Particle::get_mass() {

    return get_n()* *s->findNearest<Mass>()->get_data();
}


double Particle::get_volume()  {

    return get_n() * nt->get_m_volume();
}


double Particle::get_diameter() {

    return 1.240700981798800 * pow(get_volume(), 1.0/3.0);
}


double Particle::get_surface() {

    return 4.835975862049408 * pow(get_volume(), 2.0/3.0);
}

#endif // PARTICLE_H
