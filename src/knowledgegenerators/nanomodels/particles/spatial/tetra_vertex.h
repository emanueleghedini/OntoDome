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

#ifndef VERTEX_H
#define VERTEX_H

#include <valarray>
#include <list>
#include <climits>

#include "../dynamic/dynamicparticle.h"

class Face; // forward declaration
class Tetrahedron; // forward declaration

class Vertex
{
    /// vertex local index
    int local_id;

    /// vertex global index
    int global_id;

    /// vertex position in the space
    std::valarray<double>* x;

    /// Pointer to the particle
    std::shared_ptr< DynamicParticle > p;

    /// vertex old position in the space (constraints algorithms)
    std::valarray<double>* x_old;


    /// List of tetrahedra faces in which the vertex is included
    std::list<Face*> faces;


    /*--------- Graph visit data --------*/

    /// status for the visit
    bool visited;

    /// Visit predecessor
    Vertex* pred;

    /// hops from the source
    int hops;

    /// Base tetrahedron
    bool base;

    // Added to the basic tetrahedron
    bool added;

public:

    /// Blank Constructor
    Vertex();

    /*/// Contructor, takes vertex position
    /// \param _x : spatial coordinates for the vertex
    /// \param _global_index: vertex's global index in the simulation
    Vertex(std::valarray<double> _x, int _global_index, int _local_index);*/

    /// Contructor, takes vertex position
    /// \param _x : spatial coordinates for the vertex
    /// \param _global_index: vertex's global index in the simulation
    Vertex(std::shared_ptr< DynamicParticle >_p, int local_index);


    /*/// Contructor, takes vertex position
    /// \param _x : spatial coordinates for the vertex
    Vertex(std::valarray<double> _x);*/

    /// Returns vertex position (reference)
    const std::valarray<double> get_x() { return *x; }

    /// Returns vertex position at previous update (reference)
    const std::valarray<double> get_x_old() { return *x_old; }

    /// Returns vertex's local id
    int get_local_id() const { return local_id; }

    /// Returns vertex'a global id
    int get_global_id() const { return global_id; }

    /// Set the flag that the vertice is added to the rigid body
    /// \param _ad: true if the vertex is in the constrained body or not
    void set_added(bool _ad) { added = _ad; }

	/// Returns the mass of the particle
	double get_mass() const { return p->get_mass(); }

	/// Returns the diameter of the particle
	double get_diameter() const { return p->get_diameter(); }

    /// Set the new position
    /// \param x_new: new vertex's coordinates
    void set_x(std::valarray<double>* x_new) {

        x_old = x;
        x = x_new;
    }

    /// Update dx
    void correct_position(std::valarray<double> dx) {
                /*std::valarray<double> new_x = *x;
                new_x += dx;
                p->set_x(new_x);*/
                p->add_dx(dx);

        }

    /// Resets vertex visit data
    void reset_vertex(int new_index) {

        visited = base=added=false;
        hops = INT_MAX;
        pred = NULL;
        faces.clear();
        local_id = new_index;

    }


    /// Add face in which the vertex is included
    /// \param _f: pointer to the face in which the vertex has been included
    void add_face(Face* _f) { faces.push_back(_f); }

    /// Returns the faces in which is included the vertex
    std::list<Face*> get_faces()const { return faces; }

    /// Returns if the vertex is in the first tetrahedron or not
    bool is_base()const { return base; }

    /// Sets if the vertex is in the first tetrahedron or not
    void set_base(bool _b) { base = _b; }

    /// Returns if the vertex is in the main body or not
    bool is_added() const { return added; }

    // VISITS Methods

    /// sets vertex as visited
    /// \param _visited: true if the vertex has been visited or not
    void set_visited(bool _visited) { visited = _visited; }

    /// Returns if the vertex is visited or not
    bool is_visited() const { return visited; }

    /// Sets the predecessor in the visit
    /// \param _p: Pointer to the predecessor vertex
    void set_pred(Vertex* _p) { pred = _p; }

    /// Returns the vertex predecessor in the visit
    Vertex* get_pred() const { return pred; }

    // Adds 1 hop from the source of the visit
    void set_hops(int _step) { hops = _step; }

    /// Returns number of hops from the visit source
    int get_hop() const { return hops; }

};

#include <climits>

Vertex::Vertex() :
	x(NULL),
	x_old(NULL),
	p(NULL),
	local_id(0), global_id(0),
	visited(false), hops(INT_MAX), base(false), added(false) { }

Vertex::Vertex(std::shared_ptr< DynamicParticle >_p, int local_index):
local_id(local_index), visited(false), hops(INT_MAX),
base(false), added(false), p(_p) {

	x = &_p->x;
	x_old = &_p->x0;
	global_id = _p->get_id();

}


//Vertex::Vertex(std::valarray<double> _x, int _global_index, int _local_index) :
//	x(_x), x_old({ 0.0, 0.0, 0.0 }),
//	local_id(_local_index), global_id(_global_index),
//	visited(false), hops(INT_MAX), base(false), added(false) {}

//Vertex::Vertex(std::valarray<double> _x) :
//	x(_x), x_old({ 0.0, 0.0, 0.0 }),
//	local_id(0), global_id(-1),
//	visited(false), hops(INT_MAX), base(false), added(false) {}

#endif // !VERTEX_H
