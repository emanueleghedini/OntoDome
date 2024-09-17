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

#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "../../../ontodome.h"

/// Container Class for Nucleation Theory models
class NucleationTheory : public SoftwareModel {

public:
    NucleationTheory() {};

    std::string getClassName() const { return "Nucleation Theory Models Container"; }

    virtual double nucleation_rate() = 0;
    virtual double stable_cluster_number() = 0;
    virtual double stable_cluster_diameter() = 0;
    virtual double condensation_rate() = 0;
    virtual double get_m_volume() = 0;
    virtual double get_m_surface() = 0;

};

#endif // NUCLEATION_H
