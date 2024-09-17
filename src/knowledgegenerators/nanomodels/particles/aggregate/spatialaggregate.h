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

#ifndef SPATIALAGGREGATE_H
#define SPATIALAGGREGATE_H

#include "fractalaggregate.h"

#include <valarray>

/// An aggregate type for which the geometry is described by the particle positions.
/// A restructuring SHAKE algorithm is implemented for rearranging the particles positions
/// after a change in the bonds constraints (e.g. for sintering).
///
/// Concepts:
/// P --> PointParticle

template<typename P>
class SpatialAggregate : public FractalAggregate<P> {

public:

    /// Constructor.
    /// \param p0 pointer to the newly created particle.
    SpatialAggregate(std::shared_ptr<P> p0) : FractalAggregate<P>(p0) { }

    /// Constructor from particles and bond lists
    SpatialAggregate(std::list<std::shared_ptr<P>> _particles,
        std::list<std::shared_ptr<ParticleBond<P>>> _bonds) :
        FractalAggregate<P>(_particles, _bonds) { }

    /// Copy Constructor
    SpatialAggregate(const SpatialAggregate& _aggregate) :
        FractalAggregate<P>(_aggregate) {}

    /// Calculate the center of mass position [m]
    std::valarray<double> get_center_of_mass() const;

    /// Aggregate fractal dimension
    double get_fractal_dimension();

    /// Get the collision diameter [m]
    double get_collision_diameter() const;

    /// Get the diameter of the smallest sphere (centered in the center of mass)
    /// that encompasses completely the aggregate [m]
    double get_enclosing_sphere_diameter() const;

    /// Check collision between two aggregates
    /// !!! ALTERNATIVE: FRIEND FUNCTION
    /// \param a0 Other colliding aggregate
    /// \param c_dist maximum colliding distance
    ///	\param std::shared_ptr<P>& p0 pointer to the colliding particles in the aggregate a0
    ///	\param std::shared_ptr<P>& p1 pointer to the colliding particles in the aggregate a1
    bool check_collision(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist,
                         std::shared_ptr<P>& p0, std::shared_ptr<P>& p1);

    /// Check collision between two aggregates (periodic conditions)
    /// \param a0 Other colliding aggregate
    /// \param c_dist maximum colliding distance
    ///	\param std::valarray<double> shift coordinates shift in case of periodic conditions
    ///	\param std::shared_ptr<P>& p0 pointer to the colliding particles in the aggregate a0
    ///	\param std::shared_ptr<P>& p1 pointer to the colliding particles in the aggregate a1
    bool check_collision_periodic(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist, std::valarray<double> shift,
                                  std::shared_ptr<P>& p0, std::shared_ptr<P>& p1);

    /// Aggregate radius of gyration [m]
    double get_radius_of_gyration() const;

    /// SHAKE algorithm for bonds constraints implementation
    void shake();
};

template<typename P>
std::valarray<double> SpatialAggregate<P>::get_center_of_mass() const {

    double mass = 0;
    std::valarray<double> x = {0,0,0};

    for(auto& p: this->particles) {

        double m = p->get_mass();

        x += m*p->get_x();
        mass += m;
    }

    return x /= mass;
}


template<typename P>
double SpatialAggregate<P>::get_fractal_dimension() {

    double nr = this->get_reduced_particles_number();
    double dc = this->get_collision_diameter();
    double davg = this->get_particles_mean_diameter();

    if (dc == 0.0) // aggregate with only one particle
        return 3.0;
    else
        return log(nr)/log(dc/davg);
}


template<typename P>
double SpatialAggregate<P>::get_collision_diameter() const {

    // diameter of a sphere with the same radius of gyration of the aggregate
    return 2.0 * sqrt(5.0/3.0) * get_radius_of_gyration();
}


template<typename P>
double SpatialAggregate<P>::get_radius_of_gyration() const {

    double R_gyr = 0;
    double mass = 0;

    if (this->get_particles_number() == 1)
        return 0.0;

    std::valarray<double> x_cm = get_center_of_mass();

    for(auto& p: this->particles) {

        double m = p->get_mass();

        std::valarray<double> x = (p->get_x() - x_cm);
        double x2 = (x*x).sum();

        R_gyr += m*x2;

        mass += m;
    }

    return sqrt(R_gyr/mass);
}


template<typename P>
double SpatialAggregate<P>::get_enclosing_sphere_diameter() const {

    double diam;
    std::valarray<double> r_cm;

    std::valarray<double> x_cm = get_center_of_mass();

    double enclosing_diameter = 0;

    for(auto& p: this->particles) {

        r_cm = (p->get_x() - x_cm);
        diam = 2.*sqrt((r_cm*r_cm).sum()) + p->get_diameter();

        if(diam>enclosing_diameter)
            enclosing_diameter = diam;
    }

    return enclosing_diameter;
}


template<typename P>
bool SpatialAggregate<P>::check_collision(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist,
                                          std::shared_ptr<P>& _p0, std::shared_ptr<P>& _p1) {

    bool collision = false;

    if (this->get_id() != a0->get_id()) {

        std::valarray<double> diff = this->get_center_of_mass() - a0->get_center_of_mass();

        // Aggregate distance
        double d = sqrt((diff*diff).sum());

        // CHANGE THE DIAMETER FUNCTION
        double d1 = this->get_enclosing_sphere_diameter();
        double d2 = a0->get_enclosing_sphere_diameter();

        if (d <= c_dist + 0.5*(d1 + d2)) {
            // Check if two particles of the two aggregates are colliding properly
            for (auto p0 = this->particles.begin(); p0 != this->particles.end(); p0++) {
                for (auto p1 = a0->particles.begin(); p1 != a0->particles.end(); p1++) {
                    std::valarray<double> diff = (*p0)->get_x() - (*p1)->get_x();
                    double dist = sqrt((diff*diff).sum());
                    double p0_d = (*p0)->get_diameter();
                    double p1_d = (*p1)->get_diameter();
                    if (dist <= c_dist + 0.5*(p0_d + p1_d)) {

#ifdef VERBOSE
                        std::cout << "COLLISION!! " << std::endl;
                        std::cout << "ID(0): " << this->get_id() << " |P|: " << this->get_particles_number() << " Diameter: " << d1 << std::endl;

                        std::cout << " ID(1): " << a0->get_id() << " |P|: " << a0->get_particles_number() << " Diameter: " << d2 << std::endl;
                        std::cout << "Contact Particles: " << (*p0)->get_id() << " " << (*p1)->get_id() << std::endl
                            << " Distance(Contact Particles): " << dist
                            << " Distance(Contact Particles Recorded)" << 0.5*(p0_d + p1_d)
                            << "[m]" << std::endl;

                        std::cout << std::endl;

                        std::cout << "Particle Masses: " << (*p0)->get_mass() << ", " << (*p1)->get_mass() << std::endl;
                        std::cout << std::endl;

#endif
                        _p0 = (*p0);
                        _p1 = (*p1);
                        collision = true;
                    }
                }
            }
        }
    }

    return collision;
}
template<typename P>
bool SpatialAggregate<P>::check_collision_periodic(std::shared_ptr<SpatialAggregate<P>> a0, double c_dist,
                                                   std::valarray<double> shift,
                                                   std::shared_ptr<P>& _p0, std::shared_ptr<P>& _p1) {
    bool collision = false;

    if (this->get_id() != a0->get_id()) {

        /// get centers of mass of the aggregates
        std::valarray<double> cm1 = this->get_center_of_mass();
        std::valarray<double> cm2 = a0->get_center_of_mass();

        // shift coordinates of checked aggregate (a0)
        cm1 += shift;

        // Aggregates distance
        std::valarray<double> diff = cm1 - cm2;
        double d = sqrt((diff*diff).sum());

        double d1 = this->get_enclosing_sphere_diameter();
        double d2 = a0->get_enclosing_sphere_diameter();

        if (d <= c_dist + 0.5*(d1 + d2)) {
            // Check if two particles of the two aggregates are colliding properly
            for (auto p0 = this->particles.begin(); p0 != this->particles.end(); p0++) {
                for (auto p1 = a0->particles.begin(); p1 != a0->particles.end(); p1++) {
                    std::valarray<double> cm_p1 = (*p0)->get_x();
                    std::valarray<double> cm_p2 = (*p1)->get_x();
                    cm_p2 += shift;
                    std::valarray<double> diff = cm_p1  - cm_p2;
                    double dist = sqrt((diff*diff).sum());
                    double p0_d = (*p0)->get_diameter();
                    double p1_d = (*p1)->get_diameter();
                    if (dist <= c_dist + 0.5*(p0_d + p1_d)) {

#ifdef VERBOSE
                        std::cout << "COLLISION!! " << std::endl;
                        std::cout << "ID(0): " << this->get_id() << " |P|: " << this->get_particles_number() << " Diameter: " << d1 << std::endl;

                        std::cout << " ID(1): " << a0->get_id() << " |P|: " << a0->get_particles_number() << " Diameter: " << d2 << std::endl;
                        std::cout << "Contact Particles: " << (*p0)->get_id() << " " << (*p1)->get_id() << std::endl
                            << " Distance(Contact Particles): " << dist
                            << " Distance(Contact Particles Recorded)" << 0.5*(p0_d + p1_d)
                            << "[m]" << std::endl;

                        std::cout << std::endl;

                        std::cout << "Particle Masses: " << (*p0)->get_mass() << ", " << (*p1)->get_mass() << std::endl;
                        std::cout << std::endl;

#endif
                        _p0 = (*p0);
                        _p1 = (*p1);
                        collision = true;
                    }
                }
            }
        }
    }

    return collision;

}


template<typename P>
void SpatialAggregate<P>::shake() {

    const double SHAKE_STABILIZATION_THRESHOLD = 1e-20;

    for(auto& b: this->bonds) {

        double d = b->get_d();

        double m0 = b->get_p0()->get_mass();
        double m1 = b->get_p1()->get_mass();

        std::valarray<double> x01 = b->get_p1()->get_x() - b->get_p0()->get_x();
        std::valarray<double> x01_old = b->get_p1()->get_x_old() - b->get_p0()->get_x_old();

        double r = d*d - (x01*x01).sum();

        double den = (x01_old*x01).sum();

        // check if divide by zero
        if(fabs(den)<SHAKE_STABILIZATION_THRESHOLD) {
            double sig = (den<0) ? -1 : 1;
            den = sig*SHAKE_STABILIZATION_THRESHOLD;
        }

        double lambda = r/((1./m0 + 1./m1)*den);

        // apply the corrections to p0 and p1 positions
        b->get_p0()->add_dx(-x01_old*lambda*0.5/m0);
        b->get_p1()->add_dx(+x01_old*lambda*0.5/m1);
    }
}

#endif // SPATIALAGGREGATE_H
