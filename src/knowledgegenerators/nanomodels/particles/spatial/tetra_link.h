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

#ifndef LINK_H
#define LINK_H

#include "tetra_vertex.h"
//#include "constants.h"

#include <utility>


class Link {

protected:
	/// edge's index in the adjacency matrix
	std::pair<int, int> e_index;

	/// tuple of pointers to the edge's extemes 
	std::pair<Vertex*, Vertex*> edge;

	/// Flag for indicating if is in the first tetrahedron
	bool base;

	/// Debug variable
	double stable_lenght;

public:

	/// Default Constructor
	Link();

	/// Constructor, edge index and couple of extremes
	Link(std::pair<Vertex*, Vertex*> _edge);

	/// Returns edge extremes
	std::pair<Vertex*, Vertex*> get_edge() const { return edge; }

	/// Returns edge indices
	std::pair<int, int> get_indices() const { return e_index; }

	/// Returns the first edge's extreme
	Vertex* get_first() const { return edge.first; }

	/// Returns the second edge's extreme
	Vertex* get_second() const { return edge.second; }

	/// Sets base member
	void set_base(bool _b) { base = _b; }

	/// Sets extremes
	void set_extremes(std::pair<Vertex*, Vertex*> p) { edge = p; }

	/// Returns base member value
	bool is_base() const { return base; }

	/// Returns the lenght of the edge
	virtual double get_dist() = 0;

	/// Shake constraint evaluation
	double rearrange();

	/// Debug function 
	double get_stable_dist() { return stable_lenght; }

	/// Debug function
	void set_stable_dist(double _dist) { stable_lenght = get_dist(); }
};

#include<iostream>

//#define VERBOSE



Link::Link() {}

Link::Link(std::pair<Vertex*, Vertex*> _edge):
	edge(_edge),
	e_index(_edge.first->get_local_id(), _edge.second->get_local_id()),
	base(false){}

double Link::rearrange() {

	double D_THS = 1.0e-22;
	double Dist_THS = 1.0e-10;

	double diff = 0.0;

	double r1, r2;
	r1 = edge.first->get_diameter()*0.5;
	r2 = edge.second->get_diameter()*0.5;


	double half_diameters = r1 + r2;

	double dist = get_dist();



	//double m1, m2; m1 = m2 = 1.0; // dummy mass value
	double m1 = edge.first->get_mass();
	double m2 = edge.second->get_mass();

	std::valarray<double> x1 = edge.first->get_x() - edge.second->get_x();
	std::valarray<double> x1_old = edge.first->get_x_old() - edge.second->get_x_old();

#ifdef NANO_DEBUG
	std::valarray<double> xi = edge.first->get_x();
	std::valarray<double> xj = edge.second->get_x();
	std::valarray<double> check = edge.first->get_x() - edge.second->get_x();
	double check_dist = sqrt((check*check).sum());
	//std::cout << "R1: " << r1 << " R2: " << r2 << std::endl;
	std::cout << "Edge: " << get_indices().first << " " << get_indices().second << " "
		<< "Constraint lenght: " << dist << std::endl;
	std::cout << "ACTUAL DISTANCE: " << check_dist << std::endl;
	std::cout << "ACTUAL Xi: " << xi[0] << " " << xi[1] << " " << xi[2] << std::endl;
	std::cout << "ACTUAL Xj: " << xj[0] << " " << xj[1] << " " << xj[2] << std::endl;
	std::valarray<double> check_old = edge.first->get_x_old() - edge.second->get_x_old();
	double dist_old = sqrt((check_old*check_old).sum());
	std::cout << "OLD DISTANCE: " << dist_old << std::endl;
	std::cout << "MASSES: " << m1 << " " << m2 << std::endl;
#endif


	double r =  dist*dist - (x1*x1).sum();

	std::valarray<double> y = x1_old*x1;

	double den = (y).sum();


	if (fabs(den)<D_THS) {
		//std::cout << "SHAKE DEMONIMATOR STABILIZATION " << den << std::endl;
		double sig = (den<0) ? -1 : 1;
		den = sig*D_THS;
		//std::cout << " " << den << std::endl;
	}

	//if (fabs(r)<T_THS) {
	//	//std::cout << "SHAKE NOMINATOR STABILIZATION " << r << std::endl;
	//	r = 0.0;
	//	//std::cout << " " << den << std::endl;
	//}
	double lambda = r / ((1. / m1 + 1. / m2)*den);

	std::valarray<double> dx1, dx2;

	dx1 = +x1_old*lambda*0.5 / m1;
	dx2 = -x1_old*lambda*0.5 / m2;

	//std::valarray<double> new_x1 = edge.first->get_x() + dx1;
	//std::valarray<double> new_x2 = edge.second->get_x() + dx2;

	//edge.first->set_x(new_x1);
	//edge.second->set_x(new_x2);
	edge.first->correct_position(dx1);
	edge.second->correct_position(dx2);


	//x1 = edge.first->get_x() - edge.second->get_x();

	diff = fabs(dist*dist - (x1*x1).sum()) / (2 * dist*dist);



#ifdef NANO_DEBUG
	std::valarray<double> act_diff = edge.first->get_x() - edge.second->get_x();
	double act_dist = sqrt((act_diff*act_diff).sum());

	x1 = edge.first->get_x() - edge.second->get_x();

	std::cout<<"Edge: "<<get_indices().first<<" "<<get_indices().second<<" "
		<< "Constraint lenght: " << dist << " MODIFIED DISTANCE: " << act_dist << std::endl;
	std::cout << "X OLD: " << x1_old[0] << " " << x1_old[1] << " " << x1_old[2] << std::endl;

	std::cout << "MODIFIED Xi: " << xi[0] << " " << xi[1] << " " << xi[2] << std::endl;
	std::cout << "MODIFIED Xj: " << xj[0] << " " << xj[1] << " " << xj[2] << std::endl;

	std::cout << "MASSES: " << m1 << " " << m2 << std::endl;
	std::cout << "R: " << r << " " << "DEN: " << den << std::endl;
	std::cout << "LAMBDA: " << lambda << std::endl;
	std::cout << "DIFF: " << diff << std::endl;
	std::cout << "DX1: " << (dx1[0]*m1) << " " << (dx1[1]*m1) << " " << (dx1[2]*m1) << std::endl;
	std::cout << "DX2: " << (dx2[0]*m2) << " " << (dx2[1]*m2) << " " << (dx2[2]*m2) << std::endl;
	std::cout << std::endl;
#endif


	return diff;

}

#endif // !LINK_H
