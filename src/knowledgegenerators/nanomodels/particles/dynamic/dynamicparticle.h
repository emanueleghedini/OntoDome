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

#ifndef DYNAMICPARTICLE_H
#define DYNAMICPARTICLE_H

#include "../particle/particle.h"
#include "dynamicpoint.h"


class DynamicParticle : public Particle, public DynamicPoint {

public:

    /// Constructor
    ///	\param double _n number of molecules
    ///	\param Species _s Species Molecules Species
    ///	\param std::valarray<double> _x Molecule position
    DynamicParticle(double _n, SingleComponentComposition* _s, std::valarray<double> _x) :
        Particle(_n,_s,nt), DynamicPoint(_x) {}


};

#endif // DYNAMICPARTICLE_H
