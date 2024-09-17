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

#ifndef OBJECTCOUNTER_H
#define OBJECTCOUNTER_H


/// A class in which the id for the object is generated and used to provide a
/// common identity operator for all the derived classes.

#include <iostream>

template<typename T>
class ObjectCounter {

	static int objects_number; ///< Global counter for T object

    int id; ///< Object univocal id

public:

    /// Constructor. The counter is incremented by one.
    ObjectCounter() { 
		this->id = (this->objects_number)++;
	}

    ObjectCounter(const ObjectCounter&) { this->id = (this->objects_number)++; }

    /// Get object univocal id.
    int get_id() const { return id; }

    /// Equivalence operator.
    bool operator==(const ObjectCounter& p1) const {
        return (id==p1.id) ? true : false;
    }
};


template<typename T> int ObjectCounter<T>::objects_number(0);

#endif // OBJECTCOUNTER_H
