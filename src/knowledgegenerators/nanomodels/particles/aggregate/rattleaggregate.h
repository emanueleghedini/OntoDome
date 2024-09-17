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

#ifndef RATTLEAGGREGATE_H
#define RATTLEAGGREGATE_H

#include "spatialaggregate.h"
#include "pbmaggregate.h"
#include "../spatial/tetra_constrainer.h"
//#include "gasphase/gasphase.h"

/// Concept:
/// P --> DynamicParticle


template<typename P>
class RATTLEAggregate : public SpatialAggregate<P> {

    /*/// Pointer to the aggregate that has cannibalized the aggregate (COAGULATION)
    RATTLEAggregate* coagulated;*/

    /// Graph Support class
    Constrainer graph;

    /*/// flag if the aggregate has sintered or not
    bool sintered;*/

    /// Actual number of particles in the aggregate (sintering)
    int actual_particles;

public:

    /// Constructor
    ///	\param:	std::shared_ptr<P> p0 shared pointer to the first particle
    RATTLEAggregate(std::shared_ptr<P> p0)
        :SpatialAggregate<P>(p0), graph(p0), actual_particles(1) { }

    /// Constructor
    ///	\param:	std::list < std::shared_ptr<P> _particles list of pointers to particles
    ///	\param: std::list<std::shared_ptr<ParticleBond<P>>> _bonds list of pointers to physical bonds
    RATTLEAggregate(std::list < std::shared_ptr<P>> _particles,
        std::list<std::shared_ptr<ParticleBond<P>>> _bonds):
        SpatialAggregate<P>(_particles, _bonds) {

        // Init graph data structures
        graph.add_particles(this->particles, this->bonds);

        // Create graph
        graph.create_constraints();

        // actual particles
        actual_particles = _particles.size();
    }

    /// Copy Constructor
    ///	\param const RATTLEAggregate& _aggregate: Data to copy
    RATTLEAggregate(const RATTLEAggregate& _aggregate):SpatialAggregate<P>(_aggregate) {
        // Initialize graph
        graph.add_particles(this->particles, this->bonds);
        // Create contraints
        graph.create_constraints();
        // actual particles
        actual_particles = this->particles.size();
    }

    /// Langevin x update
    /// \param dt timestep [s]
    /// \param gamma friction coefficient [1/s]
    /// \param T temperature [K]
    void langevin_x_update(double dt, GasModels* gp, double T);

    /// Langevin v update
    /// \param dt timestep [s]
    void langevin_v_update(double dt);

    /// Initialize velocity
    ///	\param
    void v_initialization(double T);

    /// Langevin F update
    /// \param dt timestep [s]
    void langevin_F_update(double dt);

    /// Wall Reflection Management
    /// !!! POSSIBLE ALTERNATIVE: FRIEND FUNCTION
    void wall_reflection(double cube_side);

    /// Periodic Boundaries
    void periodic_domain(double _v_side_lenght);

    /// Shifts coordinates of the entire aggregate
    void shift_coordinates(std::valarray<double> _shift);

    /// Scales aggregate Coordinates
    void scale_coordinates(std::valarray<double> _scale);

    /// Change aggregate's center of mass coordinates
    void change_center_of_mass(std::valarray<double> _c_mass);

    /// Updates Constraints Graph
    void update_graph();

    /// Updates actual number of particles
    void update_particles_cardinality();

    /// Erase graph data structure
    void erase_graph();

    /// SHAKE algorithm
    bool SHAKE();

    /// Checks if some particles in the aggregate has coalesced
    bool check_coalescence();

    /// Updates constraints values
    void update_constraints();

    /// Returns the graph
    Constrainer get_graph() { return graph; }

    /// Returns the list of particles composing the aggregate
    /// !!! CHANGE APPROACH: PASS TO THE CLASS A CONTAINER TO FILL WITH DESIRED DATA OR A FUNCTION FOR FEEDING SOMETHING.
    std::list<std::shared_ptr<P>> get_particles() const;

    /// Test Function
    void print_test();

};

template<typename P>
void RATTLEAggregate<P>::langevin_x_update(double dt, GasModels* gp, double T) {

    // for all DynamicParticles composing the aggregate, update x
    for (auto part: this->particles) {

        double dp = part->get_diameter();

        // free molecular friction coefficient [1/s] (in the reference article alpha = mass*gamma)
        double gamma = (4./3.) * M_PI * gp->get_gas_flux() * dp*dp / (4.*part->get_mass());

        part->langevin_verlet_x(dt, gamma, T);
    }
}


template<typename P>
void RATTLEAggregate<P>::langevin_v_update(double dt) {

    // for all DynamicParticles composing the aggregate, update v
    for (auto part : this->particles) {
        part->langevin_verlet_v(dt);
    }

}

template<typename P>
void RATTLEAggregate<P>::v_initialization(double T) {

    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {
        (*p)->init_maxwellian_v(T);
    }

}

template<typename P>
void RATTLEAggregate<P>::langevin_F_update(double dt) {
    // TODO
}


template<typename P>
void RATTLEAggregate<P>::wall_reflection(double cube_side) {

    for (auto part : this->particles)
        part->wall_reflection(cube_side);
}

template<typename P>
void RATTLEAggregate<P>::periodic_domain(double _v_side_lenght) {

    for (auto part : this->particles)
        part->periodic_domain(_v_side_lenght);

}

template<typename P>
void RATTLEAggregate<P>::shift_coordinates(std::valarray<double> _shift) {

    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

        std::valarray<double> x = (*p)->get_x();
        x += _shift;
        (*p)->set_x(x);
    }

}

template<typename P>
void RATTLEAggregate<P>::change_center_of_mass(std::valarray<double> n_c_mass) {

    std::valarray<double> old_com = this->get_center_of_mass();
    std::valarray<double> diff = n_c_mass - old_com;

    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

        std::valarray<double> x = (*p)->get_x();
        x = x + diff;
        (*p)->set_x(x);
    }

}

template<typename P>
void RATTLEAggregate<P>::scale_coordinates(std::valarray<double> _scale) {

    for (auto p = this->particles.begin(); p != this->particles.end(); p++) {

        std::valarray<double> x = (*p)->get_x();
        x *= _scale;
        (*p)->set_x(x);
    }

}

template<typename P>
void RATTLEAggregate<P>::update_graph() {

#ifdef VERBOSE
    std::cout << "UPDATING GRAPH STRUCTURE" << std::endl;
#endif

    // Reset graph data
    graph.reset_graph();

    // Create new data structures
    graph.add_particles(this->particles, this->bonds);

    // Compute constraints
    graph.create_constraints();
}


template<typename P>
void RATTLEAggregate<P>::update_particles_cardinality() {
    actual_particles = this->particles.size();
}


template<typename P>
void RATTLEAggregate<P>::erase_graph() { graph.reset_graph(); }


template<typename P>
bool RATTLEAggregate<P>::SHAKE() {

    if (graph.SHAKE()) {
        return true;
    }
    else {
        std::cout << "SHAKE for Aggregate: " << this->get_id()<<" Diverged, ELIMINATED!!"<<std::endl;
        return false;
    }
}


template<typename P>
bool RATTLEAggregate<P>::check_coalescence() {
    if(actual_particles > this->particles.size() && this->particles.size()>0) {

#ifdef VERBOSE
        std::cout << "AGGREGATE ID: " << this->get_id() << std::endl;
        std::cout << "Actual Particles: " << actual_particles << std::endl;
        std::cout << "Listed Particles: " << this->particles.size() << std::endl;
#endif

        return true;

    } else {

        return false;
    }
}


template<typename P>
void RATTLEAggregate<P>::update_constraints() {
    graph.update_bonds(this->bonds);
}


template<typename P>
std::list<std::shared_ptr<P>> RATTLEAggregate<P>::get_particles()const { return this->particles; }


template<typename P>
void RATTLEAggregate<P>::print_test() {

    // print particle position
    for (auto it = this->particles.begin(); it != this->particles.end(); ++it) {
        (*it)->print_point();
    }
    // print vertices position of the graph structure
    //this->graph.print();
}

#endif // RATTLEAGGREGATE_H
