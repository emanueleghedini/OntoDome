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

#ifndef ONTODOME_BASE
#define ONTODOME_BASE

#include <list>
#include <vector>
#include <string>

// Useful constants
#define K_BOL 1.380650524e-23    // [J/K]
#define E_CHARGE 1.602176565e-19 // [C]
#define AMU 1.660538921e-27      // [kg]
#define N_AVO 6.0221408571e23    // [#/mol]


// Normalization units
#define MASS_UNIT 1.66053892e-027   //kg (= 1 amu)
#define LENGTH_UNIT 1.0e-010        //m (= 1 angstrom)
#define TIME_UNIT 1.01805055e-014   //s
#define ENERGY_UNIT 1.60217657e-019 //J (= 1 ev)
#define FORCE_UNIT 1.60217657e-009  //N
#define TEMPERATURE_UNIT (1.0/8.617332478e-5) //K


/// Useful math macro
#define square(x) ((x)*(x))


// activate event log
//#define VERBOSE

// activate DEBUG prints
//#define NANO_DEBUG

// Parameters used for SHAKE algorithm
#define SHAKE_CONV_CRITERION 1e-4
#define SHAKE_MAX_ITERATIONS 500

// struct for raw data from a configuration file;
typedef struct raw_c_data {

	// printing and savings/reaching data
	int SAVE_STEPS = -1;
	int PRINT_STEPS = -1;;
	int SAVE_VTK_STEPS = -1;
	std::string streams_path = "";
	std::string save_path = "";
	std::string vtk_save = "";
	std::string stream_save_path = "";
	std::list<std::string> params_to_print;
	// conditions
	std::list<std::string> c_species;
	std::list<double> c_s_molar_fraction;
	std::list<std::string> b_species;
	std::list<double> b_s_molar_fraction;
	double pressure = -1;
	double start_temp = -1;
	double end_temp = -1;
	double end_time = -1;
	double t_grad = -1;
	double volume = -1;
	//models
	std::string model_type;
	int max_part = -1;
	int min_part = -1;
	double max_time_step = -1;
	double frac_dim = -1;
	double dt = -1;
}raw_configuration_data;

#endif // ONTODOME_BASE
