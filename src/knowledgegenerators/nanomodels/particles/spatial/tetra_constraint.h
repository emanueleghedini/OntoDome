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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include "tetra_link.h"

/// Class for dummy bonds among particles to maintain the rigid body motion by means of a constrained algorithm
class Constraint : public Link
{
    /// Cosine of the angle subtended to the Contraint
    double cos_angle;

    // Subtended Vertex
    Vertex* angle_v;

    // Pointers to other edges composing the angle
    Link* p0;
    Link* p1;

public:

    /// Blank Constructor
    Constraint();

    /// Constructor
    /// \param _edge:	pair of pointers to the edge's extremes
    /// \param _p0:		pointer to the first side of the angle subtended by the edge
    /// \param _p1:		pointer to the second side of the angle subtended by the edge
    /// \param _angle:	pointer to the vertex of the angle subtended by the edge
    Constraint(std::pair<Vertex*, Vertex*> _edge, Link* _p0, Link* _p1, Vertex* _angle);

    /// Constructor
    //Constraint(std::pair<Vertex*, Vertex*> _edge);

    /// Returns Cosine of the angle
    double get_angle() const { return cos_angle; }

    /// Sets the angle subtebded by the constraint
    void set_angle(Vertex* _v, Link* _p0, Link* _p1);

    /// Returns the distance
    double get_dist();
};

#include <cfloat>


Constraint::Constraint() :
	p0(NULL), p1(NULL), angle_v(NULL), cos_angle(DBL_MAX) {
	edge = std::pair<Vertex*, Vertex*>(new Vertex(), new Vertex);
	e_index = std::pair<int, int>(-1, -1);
	base = false;
}

//Constraint::Constraint(std::pair<Vertex*, Vertex*> _edge) :
//	p0(NULL), p1(NULL), angle_v(NULL), cos_angle(DBL_MAX) {
//
//	edge = _edge;
//	e_index = std::pair<int, int>(edge.first->get_local_id(), edge.second->get_local_id());
//	base = false;
//}


Constraint::Constraint(std::pair<Vertex*, Vertex*> _edge, Link* _p0, Link* _p1, Vertex* _angle) :
	p0(_p0), p1(_p1), angle_v(_angle) {

	edge = _edge; e_index = std::pair<int, int>(edge.first->get_local_id(), edge.second->get_local_id());
	base = false;

	// Compute the cosine
	std::valarray<double> diff1 = _p0->get_first()->get_x() - _p0->get_second()->get_x();
	std::valarray<double> diff2 = _p1->get_first()->get_x() - _p1->get_second()->get_x();

	double d1 = sqrt((diff1*diff1).sum());
	double d2 = sqrt((diff2*diff2).sum());

	//cos_angle = (diff1*diff2).sum() / (d1 * d2);



	// Check distance
	std::valarray<double> check_diff = get_first()->get_x() - get_second()->get_x();
	double check_dist = sqrt((check_diff*check_diff).sum());

	//cos_angle = ((check_diff*check_diff).sum() - (diff1*diff1).sum() - (diff2*diff2).sum()) / (-2.0*d1*d2);
	cos_angle = ((diff1*diff1).sum() + (diff2*diff2).sum() - (check_diff*check_diff).sum())/ (2.0*d1*d2);

	//double cosine_dist = sqrt(d1*d1 + d2*d2 - 2.0*d1*d2*cos_angle);

}

void Constraint::set_angle(Vertex* _v, Link* _p0, Link* _p1) {

	p0 = _p0; p1 = _p1; angle_v = _v;

	// Compute the cosine
	std::valarray<double> diff1 = _p0->get_first()->get_x() - _p0->get_second()->get_x();
	std::valarray<double> diff2 = _p1->get_first()->get_x() - _p1->get_second()->get_x();

	cos_angle = (diff1*diff2).sum() / (sqrt((diff1*diff1).sum()) * sqrt((diff2*diff2).sum()));
}

double Constraint::get_dist() {

	double dist = 0.0;

	double d1 = p0->get_dist();
	double d2 = p1->get_dist();

	dist = sqrt(d1*d1 + d2*d2 - 2.0*d1*d2*cos_angle);

	return dist;


}

#endif // CONSTRAINT_H
