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

#ifndef MOMENTS_H
#define MOMENTS_H

#include "../../../ontodome.h"

/// Container Class for Moment models
class MomentsModel : public SoftwareModel {

public:
    MomentsModel() {};

    std::string getClassName() const { return "Moment Models Container"; }

    virtual double get_n_density() = 0;
    virtual double get_mean_diameter() = 0;
    virtual double get_total_volume() = 0;
    virtual double get_total_area() = 0;
    virtual double timestep(double dt) = 0;
    virtual double geometric_mean_v() = 0;
    virtual double ln2_standard_dev() = 0;
    virtual double standard_dev() = 0;
    virtual double get_M0() = 0;
    virtual double get_M1() = 0;
    virtual double get_M2() = 0;
    virtual void print_lognormal_val(const std::string path) = 0;
    virtual void get_lognormal_val(int _s_index) = 0;

};

#endif // MOMENTS_H
