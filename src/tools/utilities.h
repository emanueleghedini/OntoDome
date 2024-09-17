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

#ifndef _UTILITIES_H_
#define _UTILITIES_H_

#include <cmath>

#include <iomanip>
#include <sstream>
#include <tuple>
#include <unordered_set>
#include <memory>
#include <list>
//#include "aggregate/rattleaggregate.h"
//#include "dynamicparticle.h"

/// Simple function for triggering an output based on iterations number
bool counter_trigger(std::size_t iter, std::size_t every_n_iter);

//template<class dataType>std::string V_To_S(std::vector<dataType> vect) {


//	std::string ret_str;

//	for (auto it : vect) {

//		std::stringstream ss;

//		ss << std::setprecision(4) << it;

//		std::string part;

//		ss >> part;
//		ret_str += part + " ";
//	}

//	return ret_str;

//}

///// Functions for creating the VTK file in langevin simulations
//std::vector<double> getParticleCoords(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
//std::vector<double> getParticleMass(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
//std::vector<int>	getParticleAggregate(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
//std::vector<double> getParticleDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
//std::vector<double>	getAggregatesDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
//std::vector<int>	getParticleID(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

//int getTotalEdges(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

//std::vector<std::tuple<int, int, int>> getConnections(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

//// vect[0] -> couples
//// vect[1] -> types
//// vect[2] -> offsets
//std::vector<std::vector<int>> GetAbsoluteIndices(std::vector<std::tuple<int, int, int>> graph_data, std::vector<int> p_list);

#endif // _UTILITIES_H_
