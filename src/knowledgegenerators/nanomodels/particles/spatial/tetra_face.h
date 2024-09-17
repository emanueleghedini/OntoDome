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

#ifndef FACE_H
#define FACE_H



#include "tetra_vertex.h"
#include "../bond/edge.h"
//#include "constants.h"
#include <tuple>





class Tetrahedron; // forward declaration


// Faces counter
static int N_FACES = 0;


class Face
{
	/// face's index
	int f_index;

	/// list of vertices composing the face
	std::tuple<Vertex*, Vertex*, Vertex*> face;

	/// Plan identified by the face data
	/// equation: ax + by + cz - d = 0 parameters
	double a, b, c, d;

	/// centre of the face
	std::valarray<double> centre;

	/// Tetrahedra in which the face is included
	std::list<Tetrahedron*> tetra;

public:

	/// Blank Constructor
	Face();

	/// Constructor, passing face index and list of vertices composing the face
	/// \param _face: tuple with the pointers to the vertices of the face
	/// \param _edges: list of the edges connecting the face vertices
	Face(std::tuple<Vertex*, Vertex*, Vertex*> _face);

	/// Get distance between the centre of the face and a vertex
	/// \param _v: pointer to the vertex
	double get_distance(Vertex* _v);

	/// Add tetrahedron in which the face is
	/// \param _t: pointer to the tetrahedron
	void add_tetrahedron(Tetrahedron* _t) { tetra.push_back(_t); }

	/// Returns the tuple of pointers to vertices composing the face
	std::tuple<Vertex*, Vertex*, Vertex*> get_vertices() const { return face; }

	/// Checks if the vertes is in the face
	/// \param _v:	Pointer to vertex to check
	bool is_in_face(Vertex* _v) const {
		if (std::get<0>(face)->get_local_id() == _v->get_local_id() ||
			std::get<1>(face)->get_local_id() == _v->get_local_id() ||
			std::get<2>(face)->get_local_id() == _v->get_local_id())
			return true;
		else
			return false;
	}

	/// Returns the other two vertices excluded the one in input
	std::list<Vertex*> get_opposite(Vertex* _v);

	/// Returns the Face GLOBAL ID
	int get_face_id() const { return f_index; }

};

#include<math.h>
#include <iostream>



Face::Face() {}

Face::Face(std::tuple<Vertex*, Vertex*, Vertex*> _face) :
	face(_face)
{
	f_index = N_FACES;

	N_FACES++;

	// create plan parameters
	std::valarray<double> a1, b1;
	a1 = std::get<0>(_face)->get_x() - std::get<1>(_face)->get_x();
	b1 = std::get<1>(_face)->get_x() - std::get<2>(_face)->get_x();

	// Cross product a1 x b1 to obtain a, b, c
	a = a1[1] * b1[2] - a1[2] * b1[1];
	b = a1[2] * b1[0] - a1[0] * b1[2];
	c = a1[0] * b1[1] - a1[1] * b1[0];

	// d = ax + by + cz
	std::valarray<double> f = std::get<0>(_face)->get_x();
	d = a*f[0] + b*f[1] + c*f[2];

	// centroid
	centre = (std::get<0>(_face)->get_x() + std::get<1>(_face)->get_x() + std::get<2>(_face)->get_x()) / 3.0;

	// check if centroid is in the plain
	double check = a*centre[0] + b*centre[1] + c*centre[2];
	double toll = 0.0001;

#ifdef VERBOSE
	if (std::abs(check - d) < toll) {
		std::cout << "Face " << f_index << " centre OK!" << std::endl;
	}
#endif

}

double Face::get_distance(Vertex* _v) {

	double dist = 0.0;

	std::valarray<double> diff = _v->get_x() - centre;

	dist = std::sqrt((diff * diff).sum());

	return dist;

}

std::list<Vertex*> Face::get_opposite(Vertex* _v) {

	std::list<Vertex*> opp;

	std::list<Vertex*> tmp;
	tmp.push_back(std::get<0>(face)); tmp.push_back(std::get<1>(face)); tmp.push_back(std::get<2>(face));

	bool first, second; first = second = false;

	for (auto v = tmp.begin(); v != tmp.end(); v++) {
		if ((*v)->get_local_id() != _v->get_local_id() && !first) {
			opp.push_back((*v));
			first = true;
		}
		else if ((*v)->get_local_id() != _v->get_local_id() && !second) {
			opp.push_back((*v));
			second = true;
		}
	}

	return opp;
}

#endif // !FACE_H
