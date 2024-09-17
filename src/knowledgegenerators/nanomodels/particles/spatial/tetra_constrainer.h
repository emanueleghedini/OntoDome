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

#ifndef CONSTRAINER_H
#define CONSTRAINER_H

#include "tetra_vertex.h"
#include "tetra_constraint.h"
#include "tetra_bond.h"
#include "tetra_tetrahedron.h"
#include "tetra_link.h"
#include "tetra_constants.h"

#include "../dynamic/dynamicparticle.h"
#include "../bond/particlebond.h"


#include <vector>
#include <list>
#include <valarray>
#include <string>
#include <random>
#include <map>

//class DynamicParticle;

class Constrainer
{
    /// Mapping among local and global indices
    std::map<int, int> global_to_local;

    /// Counter for vertices local index
    int counter;

    /// Vector of the vertices composing the aggregate
    std::vector<Vertex> vertices;

    /// Adjacency matrix for the connections among vertices
    std::vector<std::vector<PhysicalLink > > bond_edges;

    /// Adjacency matrix for the contraints among vertices
    std::vector<std::vector<Constraint > > constraint_edges;

    /// List of tetrahedra composing the structure
    std::list<Tetrahedron> tetrahedra;

    /*----- Private Methods -----*/

    /// Creates contraints for a 3 vertices graph
    void constraints3();

    /// Creates constraints for a 4 vertices graph (tetrahedron)
    /// \param: std::vector<Vertex*> vertices list of vertices for creating a thetrahedron
    void constraints4(std::vector<Vertex*> vertices);

    /// Creates constraints for a graph with |V| > 4 vertices
    void constraintsN();

    /// Breadth First Search to find a connected component of the graph
    /// for creating the first tetrahedron
    /// \param index: index of the starting vertex for the search
    std::vector<Vertex*> BFS(int _index);


public:

    /// Blank Constructor
    Constrainer();

    /// Constructor
    /// \param DynamicParticle _p: first aggregate particle
    Constrainer(std::shared_ptr< DynamicParticle >_p);

    /// Create constraints
    void create_constraints();

    /// Returns number of vertices
    int get_n_vert() const { return vertices.size(); }

    /// Returns the vertices
    std::vector<Vertex> get_vertices()const { return vertices; }

    /// Returns a list with the edges
    std::list<Link*> get_edges();

    /// Returns a vector with couples of global indices for each edge (constraint and physical) and the type
    /// triple i, j, type
    std::vector<std::tuple<int, int, int>> get_edges_indices();

    /// Returns all the edges by pairs of estremes (local indices)
    std::vector<std::pair<int, int>> get_edge_extremes();

    /// Returns a list with the edges types
    std::list<std::string> get_edges_types();

    /// Returns number of edges
    int get_n_edges();

    /// Reset the Graph
    void reset_graph();

    /// Shake algorithm
    bool SHAKE();

    /// Adds Particles to the graph
    void add_particles(std::list<std::shared_ptr<DynamicParticle>> p_list,
                       std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list);

    /// Updates bond distances
    void update_bonds(std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list);

};

#include <iostream>
#include <cfloat>
#include <climits>
#include <algorithm>

bool vertex_order(std::pair<int, int> a, std::pair<int, int> b) { return a.second < b.second; }

// Returns the vertices ordered by the BFS order of visit, to guide the insertion inside the aggregate
std::list<std::pair<int, int>> get_order(std::vector<Vertex> v_list) {

    // pair: first: local index, second: hops
    std::list<std::pair<int, int>> exec_seq;

    for (auto v = v_list.begin(); v != v_list.end(); v++) {
        exec_seq.push_back(std::make_pair(v->get_local_id(), v->get_hop()));
    }

    exec_seq.sort(vertex_order);
    return exec_seq;
}

Constrainer::Constrainer():counter(0){}

Constrainer::Constrainer(std::shared_ptr< DynamicParticle >_p):counter(0) {

    Vertex v(_p, counter);
    vertices.push_back(v);
    counter++;
    // Create Adjacency matrices
    // Create adjacency matrix for bonds and constraints
    for (int i = 0; i < vertices.size(); i++) {
        std::vector<PhysicalLink> tmp_b;
        std::vector<Constraint> tmp_c;
        bond_edges.push_back(tmp_b);
        constraint_edges.push_back(tmp_c);
        for (int j = 0; j < vertices.size(); j++) {
            PhysicalLink b;
            Constraint c;
            bond_edges[i].push_back(b);
            constraint_edges[i].push_back(c);
        }
    }
}

void Constrainer::create_constraints() {

    std::vector<Vertex*> v_list;

    switch (vertices.size())
    {
    case 1:
        break;
    case 2:
        break;
    case 3:
        constraints3();
        break;
    case 4:
        for (auto v = vertices.begin(); v != vertices.end(); v++) { Vertex* p_v = &(*v); v_list.push_back(p_v); }
        constraints4(v_list);
        break;
    default:
        constraintsN();
        break;
    }


}

void Constrainer::constraints3() {

    //if (bond_edges.size() == 2) {
        for (auto v1 = vertices.begin(); v1 != vertices.end(); ++v1) {
            for (auto v2 = std::next(v1, 1); v2 != vertices.end(); ++v2) {
                int i = (*v1).get_local_id();
                int j = (*v2).get_local_id();
                if(bond_edges[i][j].get_dist() == DBL_MAX ||
                    bond_edges[j][i].get_dist() == DBL_MAX){
                    // Find the angle vertex
                    Vertex* v = NULL;
                    Link *p1, *p2; p1 = p2 = NULL;
                    for (int x = 0; x < vertices.size(); x++) {
                        if (bond_edges[i][x].get_dist() != DBL_MAX &&
                            bond_edges[j][x].get_dist() != DBL_MAX) {
                            v = &vertices[x];
                            p1 = &bond_edges[i][x];
                            p2 = &bond_edges[j][x];
                            break;
                        }

                    }
                    std::pair<Vertex*, Vertex*> e_p(&(*v1), &(*v2));
                    // Create constraint
                    Constraint c(e_p, p1, p2, v);
                    // Update data structures
                    constraint_edges[c.get_indices().first][c.get_indices().second] =
                    constraint_edges[c.get_indices().second][c.get_indices().first] = c;
                }
            }

        }
    //}

        std::cout << "3 VERTICES BODY CREATION" << std::endl;

}

void Constrainer::constraints4(std::vector<Vertex*> vertices_set) {
    /*
    for (auto v1 = vertices_set.begin(); v1 != vertices_set.end(); ++v1) {
        for (auto v2 = std::next(v1, 1); v2 != vertices_set.end(); ++v2) {
            int i = (*v1)->get_local_id();
            int j = (*v2)->get_local_id();
            if ((bond_edges[i][j].get_first() == NULL || bond_edges[j][i].get_first() == NULL) &&
                (constraint_edges[i][j].get_first() == NULL || constraint_edges[j][i].get_first() == NULL)) {
                // Find the angle vertex
                Vertex* v = NULL;
                Edge *p1, *p2; p1 = p2 = NULL;
                for (auto a = vertices_set.begin(); a != vertices_set.end(); a++) {
                    int x = (*a)->get_local_id();
                    bool side1b = bond_edges[i][x].get_first() != NULL;
                    bool side1c = constraint_edges[i][x].get_first() != NULL;
                    bool side2b = bond_edges[j][x].get_first() != NULL;
                    bool side2c = constraint_edges[j][x].get_first() != NULL;
                    if ((i != x) && (j != x) &&
                        (side1b || side1c) && (side2b || side2c)) {
                        /* if ((bond_edges[i][x].get_first() != NULL &&
                            bond_edges[j][x].get_first() != NULL) ||
                            (constraint_edges[i][x].get_first() != NULL &&
                                bond_edges[j][x].get_first() != NULL) ||
                                (constraint_edges[i][x].get_first() != NULL &&
                                    constraint_edges[j][x].get_first() != NULL)) {

                        }

                        v = &vertices[x];
                        (side1b)? p1 = &bond_edges[i][x] : p1 = &constraint_edges[i][x];
                        (side2b)? p2 = &bond_edges[j][x] : p2 = &constraint_edges[j][x];

                        std::pair<Vertex*, Vertex*> e_p((*v1), (*v2));
                        // Create constraint
                        Constraint c(e_p, p1, p2, v);
                        // Update data structures
                        constraint_edges[c.get_indices().first][c.get_indices().second] =
                            constraint_edges[c.get_indices().second][c.get_indices().first] = c;
                        break;

                    }
                    else {}

                }
            }

        }
    }
    */

    // Create angles among bonds
    std::list<Link*> tmp_e;
    int n = vertices_set.size();

#ifdef VERBOSE
    /*
    std::cout << "VERTICES FOR TETRAHEDRON" << std::endl;
    for (auto a : vertices_set) {
        std::cout << "ID: " << a->get_global_id() << " ";
    }

    std::cout << "EDGES FOR TETRAHEDRON";
    for (auto e : tmp_e) {
        std::cout << "(" << e->get_edge().first->get_global_id() << "," << e->get_edge().second->get_global_id() << ")" << " ";
    }
    std::cout << std::endl;

    std::cout << "PRINT MATRICES:" << std::endl;
    std::cout << "BONDS:" << std::endl;
    for (int i = 0; i < bond_edges.size(); i++) {
        for (int j = i; j < bond_edges[i].size(); j++) {
            if (bond_edges[i][j].get_dist() != DBL_MAX) {
                std::cout << "(" << bond_edges[i][j].get_edge().first->get_global_id() << ","
                    << bond_edges[i][j].get_edge().second->get_global_id() << ")";
            }
        }
    }
    std::cout << std::endl;
    std::cout << "CONSTRAINTS:" << std::endl;
    for (int i = 0; i < constraint_edges.size(); i++) {
        for (int j = i; j <constraint_edges[i].size(); j++) {
            if (constraint_edges[i][j].get_angle() != DBL_MAX) {
                std::cout << "(" << constraint_edges[i][j].get_edge().first->get_global_id() << ","
                    << constraint_edges[i][j].get_edge().second->get_global_id() << ")";
            }
        }
    }
    std::cout << std::endl;
    */

#endif

    /*
    for (auto v1 = vertices_set.begin(); v1 != vertices_set.end(); v1++) {
        int i = (*v1)->get_local_id();
        for (auto v2 = std::next(v1, 1); v2 != vertices_set.end(); v2++) {
            int j = (*v2)->get_local_id();
            if (bond_edges[i][j].get_dist() != DBL_MAX) {
                for (auto v3 = std::next(v2, 1); v3 != vertices_set.end(); v3++) {
                    int z = (*v3)->get_local_id();
                    if (bond_edges[j][z].get_dist() != DBL_MAX && constraint_edges[i][z].get_angle() == DBL_MAX) {
                        Constraint* c = new Constraint(std::make_pair((*v1), (*v3)),
                                                       &bond_edges[i][j], &bond_edges[j][z], *v2);
                        //tmp_e.push_back(c);
                        constraint_edges[i][z] = constraint_edges[z][i] = *c;
                    }
                }
            }
        }
    }
    */
    /*
    // Create last constraint
    for (auto v1 = vertices_set.begin(); v1 != vertices_set.end(); v1++) {
        int i = (*v1)->get_local_id();
        for (auto v2 = std::next(v1, 1); v2 != vertices_set.end(); v2++) {
            int j = (*v2)->get_local_id();
            bool c1 = constraint_edges[i][j].get_angle() != DBL_MAX;
            bool b1 = bond_edges[i][j].get_dist() != DBL_MAX;
            if (c1 || b1) { // I have a constraint or a bond from i to j
                for (auto v3 = std::next(v2, 1); v3 != vertices_set.end(); v3++) {
                    int z = (*v3)->get_local_id();
                    bool b2 = bond_edges[j][z].get_dist() != DBL_MAX;
                    bool c2 = constraint_edges[j][z].get_angle() != DBL_MAX;
                    if (b2 || c2) { // I have a bond or a constraint from j to z
                        if (bond_edges[i][z].get_dist() == DBL_MAX && constraint_edges[i][z].get_angle() == DBL_MAX) { // the link is missing between v1 and v3
                            Link *p0, *p1;
                            (b1) ? p0 = &bond_edges[i][j] : p0 = &constraint_edges[i][j];
                            (b2) ? p1 = &bond_edges[j][z] : p1 = &constraint_edges[j][z];
                            Constraint* c = new Constraint(std::make_pair((*v1), (*v3)), p0, p1, *v2);
                            //tmp_e.push_back(c);
                            constraint_edges[i][z] = constraint_edges[z][i] = *c;
                        }
                    }
                }
            }
        }
    }
    */

    for (int i = 0; i < n; i++) {
        int idx_i = vertices_set[i]->get_local_id();
        for (int j = 0; j < n; j++) {
            int idx_j = vertices_set[j]->get_local_id();
            if (j != i) {
                for (int z = 0; z < n; z++) {
                    if (z != i && z != j) {
                        int idx_z = vertices_set[z]->get_local_id();
                        bool b = bond_edges[idx_i][idx_z].get_dist() == DBL_MAX;
                        bool c = constraint_edges[idx_i][idx_z].get_angle() == DBL_MAX;
                        if (b && c) {
                            Link *p0, *p1;
                            (bond_edges[idx_i][idx_j].get_dist() != DBL_MAX) ? p0 = &bond_edges[idx_i][idx_j] : p0 = &constraint_edges[idx_i][idx_j];
                            (bond_edges[idx_j][idx_z].get_dist() != DBL_MAX) ? p1 = &bond_edges[idx_j][idx_z] : p1 = &constraint_edges[idx_j][idx_z];
                            if (p0->get_indices().first != -1 && p1->get_indices().first != -1) {
                                Constraint c(std::make_pair(vertices_set[i], vertices_set[z]), p0, p1, vertices_set[j]);
                                constraint_edges[idx_i][idx_z] = constraint_edges[idx_z][idx_i] = c;
                            }
                        }
                    }
                }
            }
        }
    }


    // Create the list of edges
    for (auto v1 = vertices_set.begin(); v1 != vertices_set.end(); v1++) {
        int i = (*v1)->get_local_id();
        for (auto v2 = std::next(v1, 1); v2 != vertices_set.end(); v2++) {
            int j = (*v2)->get_local_id();
            bool b = bond_edges[i][j].get_dist() != DBL_MAX;
            bool c = constraint_edges[i][j].get_angle() != DBL_MAX;
            if (b || c) {
                (b) ? tmp_e.push_back(&bond_edges[i][j]): tmp_e.push_back(&constraint_edges[i][j]);
            }
        }
    }


#ifdef VERBOSE
    std::cout << "VERTICES FOR TETRAHEDRON" << std::endl;
    for (auto a : vertices_set) {
        std::cout << "ID: " << a->get_global_id() << " ";
    }

    std::cout << "EDGES FOR TETRAHEDRON";
    for (auto e : tmp_e) {
        std::cout << "(" << e->get_edge().first->get_global_id() << "," << e->get_edge().second->get_global_id() << ")" << " ";
    }
    std::cout << std::endl;

    std::cout << "PRINT MATRICES:" << std::endl;
    std::cout << "BONDS:" << std::endl;
    for (int i = 0; i < bond_edges.size(); i++) {
        for (int j = i; j < bond_edges[i].size(); j++) {
            if (bond_edges[i][j].get_dist() != DBL_MAX) {
                std::cout << "(" << bond_edges[i][j].get_edge().first->get_global_id() << ","
                                 << bond_edges[i][j].get_edge().second->get_global_id() << ")";
            }
        }
    }
    std::cout << std::endl;
    std::cout << "CONSTRAINTS:" << std::endl;
    for (int i = 0; i < constraint_edges.size(); i++) {
        for (int j = i; j <constraint_edges[i].size(); j++) {
            if (constraint_edges[i][j].get_angle() != DBL_MAX) {
                std::cout << "(" << constraint_edges[i][j].get_edge().first->get_global_id() << ","
                                 << constraint_edges[i][j].get_edge().second->get_global_id() << ")";
            }
        }
    }
    std::cout << std::endl;

#endif

    std::cout << "4 VERTICES BODY CREATION" << std::endl;
    Tetrahedron t(vertices_set, tmp_e);

    tetrahedra.push_back(t);


}

void Constrainer::constraintsN() {

    // Get the first base of vertices for creating the first tetrahedron
    std::vector<Vertex*> base = BFS(0);

    // Create the first tetrahedron
#ifdef VERBOSE
    std::cout << "MULTIPLE VERTICES BODY CREATION AFTER BFS" << std::endl;
#endif
    constraints4(base);

    // Add the other vertices in the graph

    // Ordering the vertices following th BFS order(hops), to add first the external vertices that are connected to the stable structure
    std::list<std::pair<int, int>> adding_order = get_order(vertices);

    //for (auto v = vertices.begin(); v != vertices.end(); v++) {
    for(auto p = adding_order.begin(); p != adding_order.end(); p++){

        Vertex* v = &vertices[p->first];
        int v_idx = v->get_local_id();
        double min_dist = DBL_MAX;
        Face* best_face = NULL;
        Vertex* best_neighbor = NULL;

        if (!v->is_added()) {
            // for each neighbor check the faces in which is included
            for (int n_v_idx = 0; n_v_idx < vertices.size(); n_v_idx++) {
                if (bond_edges[v_idx][n_v_idx].get_dist() != DBL_MAX && vertices[n_v_idx].is_added()) {
                    Vertex *n_v = &vertices[n_v_idx];
                    std::list<Face*> faces = n_v->get_faces();
#ifdef VERBOSE/*
                    for (auto tf = faces.begin(); tf != faces.end(); tf++) {
                        std::cout << (*tf)->get_face_id();
                    }*/
#endif
                    for (auto f = faces.begin(); f != faces.end(); f++) {
                        double dist = (*f)->get_distance(&(*v));
                        // Keep the best face, the best neighborg and the best distance
                        if (dist < min_dist) { min_dist = dist; best_face = (*f); best_neighbor = n_v; }
                    }

                }
            }
            // Angle vertex is the vertex neighbor to the vertex to be inserted
            // Get the other two opposite vertices
#ifdef VERBOSE
            if (best_face == NULL || best_neighbor == NULL) {
                std::cout << "NO BEST FACE OR BEST NEIGHBOUR FOUNDED" << std::endl;
                std::cout << "VERTEX TO INSERT: " << v->get_global_id()<<std::endl;
                std::cout << "Neigbourgs: ";
                for (int i = 0; i < vertices.size(); i++) {
                    if (bond_edges[v_idx][i].get_dist() != DBL_MAX) {
                        std::cout << vertices[i].get_global_id();
                        if (vertices[i].is_added()) {
                            std::cout << "IS INSIDE!" << std::endl;
                        }
                        else {
                            std::cout << "IS NOT INSIDE!" << std::endl;
                        }
                    }

                }
                std::cout << std::endl;

                std::cout << "ADDING ORDER: ";
                for (auto p = adding_order.begin(); p != adding_order.end(); p++) {
                    std::cout << "(" <<p->first<<","<<p->second<<")"<<" ";
                }
                std::cout << std::endl;
                system("PAUSE");
                exit(0);
            }
#endif

            std::list<Vertex*> opp = best_face->get_opposite(best_neighbor);


            // Temporary list of poiters to vertices for creating the tetrahedron
            std::vector<Vertex*> tetra_v;
            for (auto i = opp.begin(); i != opp.end(); i++) { tetra_v.push_back((*i)); };
            tetra_v.push_back(&(*v)); tetra_v.push_back(best_neighbor);

            // Temporary list of pointers to edges for creating the tetrahedron
            std::list<Link*> tetra_e;

            // Create the new edges (constraints)
            for (auto v_opp = opp.begin(); v_opp != opp.end(); v_opp++) {
                if (bond_edges[(*v_opp)->get_local_id()][v->get_local_id()].get_dist() == DBL_MAX &&
                    constraint_edges[(*v_opp)->get_local_id()][v->get_local_id()].get_angle() == DBL_MAX) {
                    // Create the pointers to the angle side for the constraint

                    int opp_idx = (*v_opp)->get_local_id();
                    int v_idx = v->get_local_id();
                    int best_idx = best_neighbor->get_local_id();

                    // Side 1 the edge connecting the opposite with the v best neighbor
                    Link* p1 = NULL;

                    if (bond_edges[opp_idx][best_idx].get_dist() == DBL_MAX &&
                        constraint_edges[opp_idx][best_idx].get_angle() != DBL_MAX) {

                        p1 = &constraint_edges[opp_idx][best_idx];
                    }
                    else if(bond_edges[opp_idx][best_idx].get_dist() != DBL_MAX &&
                        constraint_edges[opp_idx][best_idx].get_angle() == DBL_MAX){

                        p1 = &bond_edges[opp_idx][best_idx];
                    }

                    // Side 2 the edge connecting v with the best neighbor
                    Link* p2 = NULL;

                    if (bond_edges[v_idx][best_idx].get_dist() == DBL_MAX &&
                        constraint_edges[v_idx][best_idx].get_angle() != DBL_MAX) {

                        p2 = &constraint_edges[v_idx][best_idx];
                    }
                    else if (bond_edges[v_idx][best_idx].get_dist() != DBL_MAX &&
                        constraint_edges[v_idx][best_idx].get_angle() == DBL_MAX ) {

                        p2 = &bond_edges[v_idx][best_idx];
                    }
                    // Create the new Contraint
                    std::pair<Vertex*, Vertex*> p_e(&(*v), (*v_opp));
                    Constraint e(p_e, p1, p2, best_neighbor);

                    // Store the constraint
                    constraint_edges[(*v_opp)->get_local_id()][v->get_local_id()] =
                    constraint_edges[v->get_local_id()][(*v_opp)->get_local_id()] = e;

                }
            }

            // create the list of edges for storing the tetrahedron
            for (auto v1 = tetra_v.begin(); v1 != tetra_v.end(); v1++) {
                int i = (*v1)->get_local_id();
                for (auto v2 = std::next(v1, 1); v2 != tetra_v.end(); v2++) {
                    int j = (*v2)->get_local_id();
                    bool b = bond_edges[i][j].get_dist() != DBL_MAX;
                    bool c = constraint_edges[i][j].get_angle() != DBL_MAX;
                    if (b || c) {
                        (b) ? tetra_e.push_back(&bond_edges[i][j]) : tetra_e.push_back(&constraint_edges[i][j]);
                    }
                }
            }


            // Create the tetrahedron and add the tetrahedron to the global list
            Tetrahedron t(tetra_v, tetra_e);

            tetrahedra.push_back(t);

            // set added to true
            v->set_added(true);
#ifdef VERBOSE
            std::cout << "ADDING VERTEX TO MULTIPLE VERTICES BODY CREATION" << std::endl;
#endif
        }
    }
}

std::vector<Vertex*> Constrainer::BFS(int _start_index) {

    std::vector<Vertex*> path;

    std::list<Vertex*> s; // Search stack
    int steps = 0; // steps from the origin

    // get the element from the list
    Vertex* start = &vertices[_start_index];
    s.push_back(start);
    //set starting vertex hops = 0
    s.back()->set_hops(0);

    while (!s.empty()) {
        Vertex *out = s.back();

        if (path.size() < 4) {
            path.push_back(out);

            out->set_base(true);
            out->set_added(true);
        }

        s.pop_back();
        if (!(out->is_visited())) {
            out->set_visited(true);
            for (int i = 0; i < vertices.size(); i++) {
                if (bond_edges[out->get_local_id()][i].get_dist() != DBL_MAX) {
                    Vertex *n = &vertices[i];
                    if (n->get_hop() == INT_MAX) {
                        n->set_hops(out->get_hop() + 1);
                        n->set_pred(out);
                        s.push_back(n);
                    }
                }

            }
        }
    }

    return path;

}

std::list<Link*> Constrainer::get_edges() {

    std::list<Link*> ret;

    int N = get_n_vert();

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (bond_edges[i][j].get_dist() != DBL_MAX)
                ret.push_back(&bond_edges[i][j]);

            if (constraint_edges[i][j].get_angle() != DBL_MAX)
                ret.push_back(&constraint_edges[i][j]);
        }
    }

    return ret;

}

std::vector<std::tuple<int, int, int>> Constrainer::get_edges_indices() {

    std::vector<std::tuple<int, int, int>> ret;

    int N = get_n_vert();

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (bond_edges[i][j].get_dist() != DBL_MAX) {
                std::tuple<int,int,int> t = std::make_tuple(bond_edges[i][j].get_first()->get_global_id(),
                                                            bond_edges[i][j].get_second()->get_global_id(),
                                                            0);
                ret.push_back(t);
            }

            if (constraint_edges[i][j].get_angle() != DBL_MAX) {
                std::tuple<int, int, int> t = std::make_tuple(constraint_edges[i][j].get_first()->get_global_id(),
                                                              constraint_edges[i][j].get_second()->get_global_id(),
                                                              1);
                ret.push_back(t);
            }

        }
    }

    return ret;

}

std::list<std::string> Constrainer::get_edges_types() {

    std::list<std::string> res;

    int N = get_n_vert();

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (bond_edges[i][j].get_first() != NULL)
                res.push_back("1");

            if (constraint_edges[i][j].get_first() != NULL)
                res.push_back("0");
        }
    }

    return res;
}

std::vector<std::pair<int, int>> Constrainer::get_edge_extremes() {

    std::vector<std::pair<int, int>> indices;

    int N = vertices.size();

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            if (bond_edges[i][j].get_first() != NULL) {
                indices.push_back(std::make_pair(bond_edges[i][j].get_edge().first->get_local_id(), bond_edges[i][j].get_edge().second->get_local_id()));
            }

            if (constraint_edges[i][j].get_first() != NULL) {
                indices.push_back(std::make_pair(constraint_edges[i][j].get_edge().first->get_local_id(), constraint_edges[i][j].get_edge().second->get_local_id()));
            }
        }
    }


    return indices;

}

int Constrainer::get_n_edges() {

    int n_edges = 0;

    int N = vertices.size();

    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            bool b = bond_edges[i][j].get_dist() != DBL_MAX;
            bool c = constraint_edges[i][j].get_angle() != DBL_MAX;
            if (b)
                n_edges++;
            if (c)
                n_edges++;

        }
    }

    return n_edges;
}

bool Constrainer::SHAKE() {

    const int SHAKE_ITERATIONS = 100;
    const double SHAKE_THS = 1.0e-4;

    int iter = 0;
    bool complete = false;

    int N = vertices.size();

    double EPS = 0.0;

    while ((iter < SHAKE_ITERATIONS) && !(complete)) {

        complete = true;

        // All the bonds
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                //if (bond_edges[i][j].get_first() != NULL) {
                if (bond_edges[i][j].get_dist() != DBL_MAX) {
                    double ths = bond_edges[i][j].rearrange();
                    bond_edges[j][i] = bond_edges[i][j];
                    EPS += ths;
                    if (ths > SHAKE_THS) {
                        complete = false;
                    }
                }
            }
        }

            // All the constraints
            for (int i = 0; i < N; i++) {
                for (int j = i + 1; j < N; j++) {
                    //if (constraint_edges[i][j].get_first() != NULL) {
                    if (constraint_edges[i][j].get_angle() != DBL_MAX) {
                        double ths = constraint_edges[i][j].rearrange();
                        constraint_edges[j][i] = constraint_edges[i][j];
                        EPS+= ths;
                        if (ths > SHAKE_THS) {
                            complete = false;
                        }
                    }
                }
            }

            iter++;
    }

    if (EPS < 1000.0) {
        return true;
    }
    else {
        std::cout << "SHAKE Epsilon: " << EPS << std::endl;
        return false;
    }

}


void Constrainer::reset_graph() {

    // Erase lists
    for (auto i = 0; i < bond_edges.size(); i++) {
        bond_edges[i].clear();
        constraint_edges[i].clear();
    }
    bond_edges.clear();
    constraint_edges.clear();

    // Erase Vertices
    vertices.clear();

    // reset counter
    counter = 0;
}

void Constrainer::add_particles(std::list<std::shared_ptr<DynamicParticle>> p_list,
    std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list) {

    //TEST
    std::cout << "particles before " << p_list.size() << std::endl;

    // Create vertices list
    for (auto part = p_list.begin(); part != p_list.end(); part++) {
        Vertex v((*part), counter);
        global_to_local[(*part)->get_id()] = counter;
        vertices.push_back(v);
        counter++;
    }

    //TEST
    std::cout << "particles after " << p_list.size() << std::endl;

    // Create Adjacency matrices
        // Create adjacency matrix for bonds and constraints
    for (int i = 0; i < vertices.size(); i++) {
        std::vector<PhysicalLink> tmp_b;
        std::vector<Constraint> tmp_c;
        bond_edges.push_back(tmp_b);
        constraint_edges.push_back(tmp_c);
        for (int j = 0; j < vertices.size(); j++) {
            PhysicalLink b;
            Constraint c;
            bond_edges[i].push_back(b);
            constraint_edges[i].push_back(c);
        }
    }
    //TEST
    std::cout << "bonds before" << b_list.size()<<std::endl;

    int a, b;
    for (auto e = b_list.begin(); e != b_list.end(); e++) {

        // get local indices
        a = global_to_local[(*e)->get_v0()->get_id()];
        b = global_to_local[(*e)->get_v1()->get_id()];

        std::shared_ptr<ParticleBond<DynamicParticle>> pb = (*e);

        // Create physical link
        PhysicalLink link(std::make_pair(&vertices[a], &vertices[b]), pb);


        bond_edges[a][b] = bond_edges[b][a] = link;
    }

    //TEST
    std::cout << "bonds after" << b_list.size() << std::endl;

}

void Constrainer::update_bonds(std::list < std::shared_ptr<ParticleBond<DynamicParticle>>> b_list) {

    for (auto e = b_list.begin(); e != b_list.end(); e++) {

        // get local indices
        int a = global_to_local[(*e)->get_v0()->get_id()];
        int b = global_to_local[(*e)->get_v1()->get_id()];

        double new_distance = (*e)->get_bond_distance();

        bond_edges[a][b].update_dist(new_distance);
        bond_edges[b][a].update_dist(new_distance);
    }

}

#endif //CONSTRAINER_H
