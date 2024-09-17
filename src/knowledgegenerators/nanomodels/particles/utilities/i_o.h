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

#ifndef I_O
#define I_O

#include <iostream>
#include <string>


class i_o{

	/// errors file path
	std::string err_file_path = "err/errors.dat";

	/// log file path
	std::string log_file_path = "log/log.dat";

public:

	/// Contructor
	i_o();

	/// function for writing to the log file
	/// \param std::string entry: entry to write
	void log_entry(std::string _entry);

	/// function for writing to the error output (blocking)
	/// \param std::string error: error description
	void error_entry_not_blocking(std::string _error);

	/// function for writing to the error output (blocking)
	/// \param std::string error: error description
	void error_entry_blocking(std::string _error);


};

i_o::i_o(){}


void i_o::log_entry(std::string entry){



}


void i_o::error_entry_not_blocking(std::string _error){

	std::cout << "ERROR!! Check errors.dat\n";
	freopen(err_file_path.c_str(), "a", stderr);
	std::cerr << _error << std::endl;
	fclose(stderr);
}

void i_o::error_entry_blocking(std::string _error){

	std::cout << "ERROR!! Check errors.dat\n";
	freopen(err_file_path.c_str(), "a", stderr);
	std::cerr << _error << std::endl;
	fclose(stderr);
	system("PAUSE");
	exit(0);

}

#endif // I_O
