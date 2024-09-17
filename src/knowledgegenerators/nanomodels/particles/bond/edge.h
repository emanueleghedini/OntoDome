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

#ifndef EDGE_H
#define EDGE_H

#include <memory>

/// This class provides the basic implementation for the edge.
template<typename V>
class Edge {

    std::shared_ptr<V> v0; ///< First vertex
    std::shared_ptr<V> v1; ///< Second vertex

public:

    /// Constructor
    /// \param _v0 the first vertex
    /// \param _v1 the second vertex
    Edge(std::shared_ptr<V> _v0, std::shared_ptr<V> _v1) : v0(_v0), v1(_v1) { }

    /// Get pointer to v0
    std::shared_ptr<V> get_v0() const { return v0; }

    /// Get pointer to v1
    std::shared_ptr<V> get_v1() const { return v1; }

    /// Set pointer to v0
    void set_v0(std::shared_ptr<V> _v0) { v0 = _v0; }

    /// Set pointer to v1
    void set_v1(std::shared_ptr<V> _v1) { v1 = _v1; }

    /// Check if an object is part of this edge
    /// \param v2 object to be checked
    bool is_vertex(std::shared_ptr<V> v2) const {
        return ((v0==v2)||(v1==v2)) ? true : false;
    }

    /// Check if an edge is adjacent
    /// \param e0 edge to be checked
    bool is_adjacent(Edge& e0) const {
        return (v0==e0->v0)||(v0==e0->v1)||
               (v1==e0->v0)||(v1==e0->v1)
                ? true : false;
    }

    /// Check if an edge is parallel
    /// \param e0 edge to be checked
    bool is_parallel(Edge& e0) const {
        return ((v0==e0->v0)&&(v1==e0->v1))||
               ((v0==e0->v1)&&(v1==e0->v0))
                ? true : false;
    }
};

#endif // EDGE_H
