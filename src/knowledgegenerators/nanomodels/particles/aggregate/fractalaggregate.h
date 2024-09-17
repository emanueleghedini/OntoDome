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

#ifndef FRACTALAGGREGATE_H
#define FRACTALAGGREGATE_H

#include "aggregate.h"
#include "../base/collisionalobject.h"


/// An aggregate type in which the geometrical complexity can be expressed in terms of
/// a fractal dimension.
///
/// Concepts:
/// P --> Particle

template<typename P>
class FractalAggregate : public Aggregate<P>, public CollisionalObject {

public:

    /// Constructor.
    /// \param p0 pointer to the newly created particle.
    FractalAggregate(std::shared_ptr<P> p0) : Aggregate<P>(p0) {}

    /// Constructor from particles and bonds lists
    ///	\param
    ///	\param
    FractalAggregate(std::list<std::shared_ptr<P>> _particles,
            std::list<std::shared_ptr<ParticleBond<P>>> _bonds) :
            Aggregate<P>(_particles, _bonds) {}

    /// Copy Constructor
    FractalAggregate(const FractalAggregate& _aggregate):Aggregate<P>(_aggregate) {}

    /// Aggregate fractal dimension
    virtual double get_fractal_dimension() = 0;
};

#endif // FRACTALAGGREGATE_H
