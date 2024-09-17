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

#ifndef TETRABOND_H
#define TETRABOND_H

#include "tetra_link.h"
#include "tetra_constants.h"
#include "../bond/particlebond.h"
#include "../dynamic/dynamicparticle.h"

#include <cfloat>

/// Class describing the contact or sintering bond between particles


class PhysicalLink : public Link {

    /* SHAKE DATA */
    /// Bonds's Lenght
    double lenght;
    //double *lenght;
    std::shared_ptr<ParticleBond<DynamicParticle>> b;

public:

    /// Blank Constructor
    PhysicalLink();

    /// Constructor
    /// \param _edge:		pair of pointers to the edge's extremes
    /// \param _bond_dist:	bond fixed distance
    //PhysicalLink(std::pair<Vertex*, Vertex*> _edge, double _bond_dist);
    PhysicalLink(std::pair<Vertex*, Vertex*> _edge, std::shared_ptr<ParticleBond<DynamicParticle>> _b);

    /// Returns the bond distance
    double get_dist() {
                if (b != NULL) {
                        lenght = b->get_bond_distance();
#ifdef NANO_DEBUG
                        std::cout << "SINTERING LEVEL: " << b->get_sintering_level() << std::endl;
#endif
                }
        else
            lenght = DBL_MAX;
                /*if (b != NULL)
                        return lenght;
                else
                        lenght = DBL_MAX;*/

        return lenght;
    }

    /// Updates bond lenght
    void update_dist(double _dist) { lenght = _dist; }

};

#include <cfloat>


PhysicalLink::PhysicalLink():
    lenght(DBL_MAX), b(NULL)
{
    e_index = std::pair<int, int>(-1, -1);
    edge = std::pair<Vertex*, Vertex*>(new Vertex(), new Vertex());
    base = false;
}

PhysicalLink::PhysicalLink(std::pair<Vertex*, Vertex*> _edge, std::shared_ptr<ParticleBond<DynamicParticle>> _b) :
    lenght(_b->get_bond_distance())
{
        b = _b;
    edge = _edge;
    e_index = std::pair<int, int>(edge.first->get_local_id(), edge.second->get_local_id());
    base = false;
}

#endif // TETRABOND_H
