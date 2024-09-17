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

#ifndef TETRAHEDRON_H
#define TETRAHEDRON_H

#include "tetra_face.h"
#include "tetra_link.h"
//#include "constants.h"

#include <list>
#include <vector>




	static int N_TETRA = 0;


	class Tetrahedron
	{
		// absolute index of the tetrahedron
		int tetra_index;

		// list of faces composing the tetrahedron
		std::list<Face> faces;

		// list of edges
		std::list<Link*> edges;

		// list of vertices
		std::list<Vertex*> vertices;


	public:

		/// Blank constructor
		Tetrahedron();

		// constructor with a set of vertex (all virtual edges)
		Tetrahedron(std::list<Vertex*> _vs);

		// constructor with a vertex and a face
		Tetrahedron(Face _f, Vertex _v);

		// constructor with a set of vertices and a set of edges
		Tetrahedron(std::vector<Vertex*> _vs,
					std::list<Link*> _e);

		// Method for getting the Edges composign the Tetrahedron
		std::list<Link*> get_edges() const { return edges; }

		// Method for getting the Vertices composing the tetrahedron
		std::list<Vertex*> get_vertices() const { return vertices; }

		// Method for getting the faces
		const std::list<Face>& get_faces() const { return faces; }

		// Get tetrahedron index
		int get_index() const { return tetra_index; }



	};

#include <iostream>
#include <utility>


Tetrahedron::Tetrahedron() {}

Tetrahedron::Tetrahedron(std::list<Vertex*> _vs) {

	//if (_vs.size() < 4 || _vs.size() > 4) {
	//	std::cout << "Not possible to create a Tetrahedron, not correct number of vertices" << std::endl;
	//}
	//else {
	//	// Create edges
	//	for (auto v1 : _vs) {
	//		for (auto v2 : _vs) {
	//			if (v1->get_id() != v2->get_id()) {
	//				std::pair<Vertex*, Vertex*> tmp(v1, v2);
	//				Edge *e = new Edge(tmp, false);
	//				edges.push_back(e);
	//			}

	//		}
	//	}

	//	// Create vertices
	//	vertices = _vs;

	//	// Create Faces
	//	for (auto v1 = _vs.begin(); v1 != _vs.end(); ++v1) {
	//		for (auto v2 = v1; v2 != _vs.end(); ++v2) {
	//			for (auto v3 = v2; v3 != _vs.end(); ++v3) {
	//				std::tuple<Vertex*, Vertex*, Vertex*> ins(*v1, *v2, *v3);
	//				faces.push_back(ins);
	//			}
	//		}
	//	}

	//	tetra_index = N_TETRA;

	//	N_TETRA++;

	//}



}

Tetrahedron::Tetrahedron(Face _f, Vertex _v) {

}

Tetrahedron::Tetrahedron(std::vector<Vertex*> _vs,
	std::list<Link*> _e) {

	if (_vs.size() == 4 && _e.size() == 6) {

		// Edges
		edges = _e;

		// vertices
		for (int i = 0; i < _vs.size(); i++) { vertices.push_back(_vs[i]); }

		// Assign index
		tetra_index = N_TETRA;
		N_TETRA++;

		// Create faces
		for (auto v1 = _vs.begin(); v1 != _vs.end(); v1++) {
			for (auto v2 = std::next(v1, 1); v2 != _vs.end(); v2++) {
				for (auto v3 = std::next(v2, 1); v3 != _vs.end(); v3++) {
					std::tuple<Vertex*, Vertex*, Vertex*> t(*v1, *v2, *v3);
					Face* f = new Face(t);
					faces.push_back(*f);

					// Update vertices informations
					(*v1)->add_face(f); (*v2)->add_face(f); (*v3)->add_face(f);
					// add face to the tetrahedron
					f->add_tetrahedron(this);
				}
			}
		}
	}
	else {
		std::cout << "NOT ENOUGH EDGES(" << _e.size() << ")" << "OR VERTICES(" << _vs.size() << ")"
			<< "FOR CREATING A TETRAHEDRON" << std::endl;
		system("PAUSE");
		exit(0);

	}

}

#endif // TETRAHEDRON_H
