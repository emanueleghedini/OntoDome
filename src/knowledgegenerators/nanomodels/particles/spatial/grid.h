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

#ifndef GRID_H
#define GRID_H

#include "../spatial/cell.h"

//#include "aggregate/spatialaggregate.h"
#include "../aggregate/rattleaggregate.h"
#include "../dynamic/dynamicparticle.h"
#include "../utilities/i_o.h"

#include <vector>


class Grid {


    /// Linked-Cell grid
    std::vector<Cell<RATTLEAggregate<DynamicParticle>>> grid;

    /// Cell side (cubic)
    double cell_side;

    /// Number of cells for side
    int C;

    /// Control Volume side
    double volume_side;

    /// Grid side
    double grid_side;

    /// Number of cells
    int n_cells;

    /// position shift from centre of the simulation(centre of the control volume) to centre of the grid (leftmost vertex)
    double shift;

    ///	Create the grid
    void grid_init();

    /// Erase the grid
    void grid_erase();

    /// Returns linear index in the grid (x, y, z)
    int get_linear(const int x1, const int x2, const int x3) {
        return x3*C*C + x2*C + x1;}

public:

    /// Blank Constructor
    Grid();

    /// Constructor
    ///	\param double _volume Control Volume dimension [m3]
    /// \param int _dimension starting number of cell for side(eg. _dimension = 5 -> grid of 5x5x5 cells)
    Grid(double _volume, int _dimension);

    /// Grid Update/Creation
    ///	\param std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble: particles ensamble
    ///	\param double _V: Control Volume dimension [m3]
    ///	\param double _agg_max_diameter: Max aggregate diameter [m]
    ///	\param int _max_id: id of the biggest particle in the ensamble [index]
    void grid_update(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble,
                     double _V, double _agg_max_diameter, int _max_id);

    /// Aggregates Dispacement on the grid
    ///	\param std::list<std::shared_ptr<SpatialAggregate<DynamicParticle>>> _ensamble: list of aggregates in the simulation
    void displace_aggregates(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble);

    /// Adds aggregates to the grid
    ///	\param std::shared_ptr<SpatialAggregate<DynamicParticle>> _p0: pointer to the aggregate to add
    void insert_aggregate(std::shared_ptr<RATTLEAggregate<DynamicParticle>> _p0 );

    /// Moves aggregates on the grid
    ///	\param std::shared_ptr<RATTLEAggregate<DynamicParticle>> _p0 pointer to the aggregate to move
    void move_aggregates();

    /// Returns total grid cells number
    int get_n_cells() const { return n_cells; }

    /// Returns grid number of cells for side
    int get_C() const { return C; }

    /// Returns the cell side of a single cell in the grid
    double get_cell_side() const { return cell_side; }

    /// Returns the grid side
    double get_grid_side()const { return grid_side; }

    /// Returns the cell of the relative coordinates
    ///	\param std::valarray<double> _x point coordinates to retrieve the cell index (_x must be a valid set of coordinates in the domain)
    std::valarray<int> get_cell_index(std::valarray<double> _x) const;

    /// Overload operator "()" for access the grid with multiple indexes (3D)
    Cell<RATTLEAggregate<DynamicParticle>> &operator()(const int _x1, const int _x2, const int _x3); //(x, y, z)

    /// Overload operator "[]" for access the grid with single index
    Cell<RATTLEAggregate<DynamicParticle>> &operator[](const int lin_index);

};

int periodic_position_grid(int c_idx, double x, int SIDE, double& shift, double domain_side) {
    if (x < -(domain_side*0.5)) { shift = +domain_side; return SIDE - 1; }
    else if (x > (domain_side*0.5)) { shift = -domain_side; return 0; }
    else {
        if (c_idx > SIDE - 1)
            return SIDE - 1;
        else if (c_idx < 0)
            return 0;
        else
            return c_idx;
    }

}

Grid::Grid():cell_side(DBL_MIN),C(0),volume_side(DBL_MIN), n_cells(0), shift(DBL_MIN), grid_side(DBL_MIN) {}

Grid::Grid(double _volume, int _dimension) {

    // Start the grid with cell side 1/10 of the control volume side
    // set the volume side
    volume_side = pow(_volume, 1.0 / 3.0);
    // update grid side
    grid_side = volume_side;
    // set the number of cells for side
    C = _dimension;
    // set the cell side
    cell_side = volume_side / (double)C;
    // set the particles shift
    shift = volume_side *0.5;
    // number of cells
    n_cells = C*C*C;

    // Initialize the grid
    grid_init();
}

void Grid::grid_update(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble,
                       double _V, double _agg_max_diameter, int _max_id) {

    // side of the control volume
    volume_side = pow(_V, 1.0 / 3.0);
    // Update shift
    shift = volume_side * 0.5;
    // difference between the control volume and the grid
    double r = abs(grid_side - volume_side);

    bool aggregate_related = cell_side < _agg_max_diameter;
    bool volume_related = r > cell_side;

    if (aggregate_related || volume_related) {
        std::cout << "Rearrange_grid" << std::endl;
        if (_agg_max_diameter > cell_side) {
            std::cout << "Aggregate Dimension" << std::endl;
            std::cout << "Aggreggate "<< _max_id <<" diameter: " << _agg_max_diameter << " Cell side: " << cell_side << std::endl;
            grid_side = volume_side;
            cell_side = 2 * _agg_max_diameter;
            C = ceil(volume_side / cell_side);
        }
        else {
            std::cout << "Volume compression" << std::endl;
            cell_side = volume_side / (double)C;
            grid_side = volume_side;
        }

        n_cells = C*C*C;
        std::cout << "Cells:" << n_cells << std::endl;

        double cell_volume = cell_side*cell_side*cell_side;
        shift = volume_side * 0.5;
        // erase the grid
        grid_erase();
        // create a new grid with new dimensions
        grid_init();
        // Displace aggregates in the grid
        displace_aggregates(_ensamble);
    }
    else {
        move_aggregates();
    }

}

void Grid::grid_erase() {

    grid.clear();
}

void Grid::grid_init() {

    for (int c = 0; c < n_cells; c++) {
        Cell<RATTLEAggregate<DynamicParticle>> cell;
        grid.push_back(cell);
    }
}

void Grid::displace_aggregates(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> _ensamble) {

    for (auto agg = _ensamble.begin(); agg != _ensamble.end(); agg++) {
        insert_aggregate((*agg));
    }

}

void Grid::insert_aggregate(std::shared_ptr<RATTLEAggregate<DynamicParticle>> _p0) {

    std::valarray<double> position = _p0->get_center_of_mass();
    double abs_x1 = position[0] + shift;
    double abs_x2 = position[1] + shift;
    double abs_x3 = position[2] + shift;

    int c_x1 = floor(abs_x1 / cell_side);
    int c_x2 = floor(abs_x2 / cell_side);
    int c_x3 = floor(abs_x3 / cell_side);

    // Periodic position displacement
    std::valarray<double> shift_coords = { 0.0,0.0,0.0 };
    c_x1 = periodic_position_grid(c_x1, position[0], C, shift_coords[0], volume_side);
    c_x2 = periodic_position_grid(c_x2, position[1], C, shift_coords[1], volume_side);
    c_x3 = periodic_position_grid(c_x3, position[2], C, shift_coords[2], volume_side);

    #ifdef VERBOSE
    if (c_x1 < 0.0 || c_x2 < 0.0 || c_x3 < 0.0) // Errors Checking
    {
        std::string error;
        error = "Negative coordinate in insert_agregate!! for point " + std::to_string(position[0]) + ", " + std::to_string(position[1]) + ", " + std::to_string(position[2]) + "\n";
        error += "SHIFT: " + std::to_string(shift) + "\n";
        error += "Shifted Coordinates " + std::to_string(position[0] + shift) + ", " + std::to_string(position[1] + shift) + ", " + std::to_string(position[2] + shift) + "\n";
        error += "Cell Coordinates " + std::to_string(c_x1) + ", " + std::to_string(c_x2) + ", " + std::to_string(c_x3) + "\n";
        error += "Id Aggregate:" + std::to_string(_p0->get_id()) + "\n";
        i_o err_grid;
        err_grid.error_entry_blocking(error);
    }
    #endif

    /// Get linear index in the grid
    int lin_index = get_linear(c_x1, c_x2, c_x3);

    /// Shift, if needed, the coordinates
    _p0->shift_coordinates(shift_coords);

    /// add the aggregate to the cell
    grid[lin_index].insert_element(_p0);
}

std::valarray<int> Grid::get_cell_index(std::valarray<double> _x) const {

    std::valarray<int> indices = { -2, -2, -2 };
    double abs_x1 = _x[0] + shift;
    double abs_x2 = _x[1] + shift;
    double abs_x3 = _x[2] + shift;

    indices[0] = floor(abs_x1 / cell_side);
    indices[1] = floor(abs_x2 / cell_side);
    indices[2] = floor(abs_x3 / cell_side);

    std::valarray<double> shift_coords = { 0.0,0.0,0.0 };

    // Periodic position displacement
    indices[0] = periodic_position_grid(indices[0], _x[0], C, shift_coords[0], volume_side);
    indices[1] = periodic_position_grid(indices[1], _x[1], C, shift_coords[1], volume_side);
    indices[2] = periodic_position_grid(indices[2], _x[2], C, shift_coords[1], volume_side);

    return indices;
}



void Grid::move_aggregates() {

    for (int idx = 0; idx < grid.size(); idx++) {

        int cell_card = grid[idx].get_n();

        for (int p_idx = 0; p_idx < cell_card; p_idx++) {


            std::shared_ptr<RATTLEAggregate<DynamicParticle>> a = grid[idx][p_idx];
            std::valarray<double> t_x = a->get_center_of_mass();

            double a_x1 = t_x[0];
            double a_x2 = t_x[1];
            double a_x3 = t_x[2];
            int c_x1 = floor((a_x1 + shift) / cell_side);
            int c_x2 = floor((a_x2 + shift) / cell_side);
            int c_x3 = floor((a_x3 + shift) / cell_side);

            std::valarray<double> shift_coords = {0.0,0.0,0.0};

            // Periodic position displacement
            c_x1 = periodic_position_grid(c_x1, a_x1, C, shift_coords[0], volume_side);
            c_x2 = periodic_position_grid(c_x2, a_x2, C, shift_coords[1], volume_side);
            c_x3 = periodic_position_grid(c_x3, a_x3, C, shift_coords[2], volume_side);

            #ifdef VERBOSE
            if (c_x1 < 0.0 || c_x2 < 0.0 || c_x3 < 0.0) // Errors Checking
            {
                std::string error;
                error = "Negative coordinate in move aggregate!! for point " + std::to_string(t_x[0]) + ", " + std::to_string(t_x[1]) + ", " + std::to_string(t_x[2]) + "\n";
                error += "SHIFT: " + std::to_string(shift) + "\n";
                error += "Shifted Coordinates " + std::to_string(t_x[0] + shift) + ", " + std::to_string(t_x[1] + shift) + ", " + std::to_string(t_x[2] + shift) + "\n";
                error += "Cell Coordinates " + std::to_string(c_x1) + ", " + std::to_string(c_x2) + ", " + std::to_string(c_x3) + "\n";
                i_o err_grid;
                err_grid.error_entry_not_blocking(error);
            }
            #endif

            int act_lin_idx = get_linear(c_x1, c_x2, c_x3);
            // The aggregate changed cell
            if (act_lin_idx != idx) {
                // Erase the aggregate for the former cell
                grid[idx].erase_n(p_idx);
                // update coodinates
                a->shift_coordinates(shift_coords);
                // insert the aggregate in the new cell
                grid[act_lin_idx].insert_element(a);
                // Update indices
                p_idx--;
                cell_card--;
            }

        }

    }

}

Cell<RATTLEAggregate<DynamicParticle>>& Grid::operator()(const int _x1, const int _x2, const int _x3) {

    int linear_idx = _x3*C*C + _x2*C + _x1;
    return grid[linear_idx];

}

Cell<RATTLEAggregate<DynamicParticle>>& Grid::operator[](const int linear_index) { return grid[linear_index]; }

#endif

