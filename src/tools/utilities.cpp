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

#include "utilities.h"


bool counter_trigger(std::size_t iter, std::size_t every_n_iter) {

    return ((std::floor(iter/every_n_iter) - double(iter)/double(every_n_iter)) == 0);
}

//// Auxiliary Function to get all the coordinates of each particle in the system
//std::vector<double> getParticleCoords(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<double> points;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {

//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			std::valarray<double> cTmp = (*p)->get_x();
//			for (int d = 0; d < 3; d++)
//				points.push_back(cTmp[d]);
//		}

//	}

//	return points;

//}

//// Auxiliary Function to get the mass of each particle in the system
//std::vector<double> getParticleMass(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<double> masses;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {

//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			masses.push_back((*p)->get_mass());
//		}
//	}

//	return masses;
//}

//// Auxiliary Function to get the ID of each nanoparticle in the system
//std::vector<int> getParticleAggregate(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	//Aggregate affiliation of particles
//	std::vector<int> aggAffil;

//	//Aggregate's affiliation id
//	int ident;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
//		ident = (*agg)->get_id();
//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			aggAffil.push_back(ident);
//		}
//	}
//	return aggAffil;


//}

//// Auxiliary Function to retrieve the ID of each particle in the aggregate
//std::vector<int> getParticleID(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<int> IDs;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
//		std::vector<int> p_IDs;
//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			p_IDs.push_back((*p)->get_id());
//		}
//		IDs.insert(IDs.end(), p_IDs.begin(), p_IDs.end());
//	}


//	return IDs;

//}

//// Auxiliary Function to get the diameter of each nanoparticle in the system
//std::vector<double> getParticleDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<double> diameter;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			diameter.push_back((*p)->get_diameter());
//		}
//	}

//	return diameter;
//}

//std::vector<double>	getAggregatesDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<double> diameter;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
//		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
//		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
//			diameter.push_back((*agg)->get_enclosing_sphere_diameter());
//		}
//	}

//	return diameter;


//}

//int getTotalEdges(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	int tot_edges = 0;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); agg++) {
//		tot_edges += (*agg)->get_graph().get_n_edges();
//	}

//	return tot_edges;

//}

//// Auxiliary function to get the connectiob data of each aggregate
//std::vector<std::tuple<int, int, int>> getConnections(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

//	std::vector<std::tuple<int, int, int>> graph_data;

//	std::vector<int> couples;
//	std::vector<int> types;
//	std::vector<int> offsets;

//	int off = 0;

//	for (auto agg = ensamble.begin(); agg != ensamble.end(); agg++) {

//		// Get graph edges (extreme indices and types)
//		std::vector<std::tuple<int, int, int>> tmp = (*agg)->get_graph().get_edges_indices();
//		graph_data.insert(graph_data.end(), tmp.begin(), tmp.end());
//	}

//	return graph_data;


//}

//// Auxiliary function to map edges extremes indices with position in the array for vtk
//std::vector<std::vector<int>> GetAbsoluteIndices(std::vector<std::tuple<int, int, int>> graph_data,
//	std::vector<int> p_list) {

//	std::vector<std::vector<int>> res;

//	std::vector<int> couples;
//	std::vector<int> types;
//	std::vector<int> offsets;

//	int off = 0;

//	for (auto t = graph_data.begin(); t != graph_data.end(); t++) {

//		//get first extreme
//		int idx1 = std::get<0>((*t));
//		int idx_abs1 = -1;
//		//get second extreme
//		int idx2 = std::get<1>((*t));
//		int idx_abs2 = -1;
//		//get type
//		int type = std::get<2>((*t));
//		// find the first index
//		for (int i = 0; i < p_list.size(); i++) {
//			if (p_list[i] == idx1) {
//				idx_abs1 = i;
//				off++;
//				break;
//			}
//		}
//		// find the second index
//		for (int i = 0; i < p_list.size(); i++) {
//			if (p_list[i] == idx2) {
//				idx_abs2 = i;
//				off++;
//				break;
//			}
//		}
//		// insert the absolute indices
//		couples.push_back(idx_abs1); couples.push_back(idx_abs2);
//		// insert the edge type
//		types.push_back(type);
//		// give offset
//		offsets.push_back(off);
//	}

//	res.push_back(couples); res.push_back(types); res.push_back(offsets);

//	return res;

//}
