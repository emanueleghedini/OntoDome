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

#include "clock.h"


WallClock::WallClock() {}


void WallClock::start() {
    t1 = std::chrono::high_resolution_clock::now();
}


void WallClock::stop() { t2 = std::chrono::high_resolution_clock::now(); }


///< evaluate interval between start and stop in seconds
double WallClock::interval() {

    std::chrono::duration<double>
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

    return time_span.count();
}

