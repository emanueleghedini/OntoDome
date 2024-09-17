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

#ifndef POINT_H
#define POINT_H

#include "../../../../ontodome.h"

#include <valarray>
#include <iostream>

/// This class provides a basic interface and data for a point object.
class Point {

    std::valarray<double> x;  ///< Object position [m]
    std::valarray<double> x0; ///< Object old position [m]

    friend class Vertex;

public:

    /// Constructor.
    /// \param _mass object mass [kg]
    /// \param _x initial position [m]
    Point(const std::valarray<double>& _x={0,0,0}) : x(_x), x0(_x) {}

    /// Get position [m].
    std::valarray<double> get_x() const { return x; }

    /// Get old position [m].
    std::valarray<double> get_x0() const { return x0; }

    /// Adjust object position.
    /// \param dx valarray containing the correction [m]
    void add_dx(const std::valarray<double>& dx) { x0 = x; x += dx; }

    /// Change object position.
    /// \param new_x new point coordinates
    void set_x(std::valarray<double> new_x) { x = new_x; }

    /// Test Function
    void print_point() {
        std::cout << "POINT:" << std::endl;
        std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
        std::cout << "POINT_OLD:" << std::endl;
        std::cout << x0[0] << " " << x0[1] << " " << x0[2] << std::endl;
    }

};

#endif // POINT_H
