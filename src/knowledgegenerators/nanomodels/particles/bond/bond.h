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

#ifndef BOND_H
#define BOND_H

#include "../../../../ontodome.h"
#include "edge.h"

#include <memory>


/// A bond is identified as a prescribed distance between two objects.
template<typename P>
class Bond : public Edge<P> {

public:

    /// Constructor
    Bond(std::shared_ptr<P> p0, std::shared_ptr<P> p1) : Edge<P>(p0,p1) {}

    /// Get the distance between particles centers [m]
    virtual double get_bond_distance() const = 0;
};

#endif // BOND_H
