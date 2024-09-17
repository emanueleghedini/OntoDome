#ifndef CELL_H
#define CELL_H


//#include <list>
#include <vector>
#include <memory>

template<typename P>
class Cell {

	/// list of aggregates in the cell
	std::vector<std::shared_ptr<P>> aggregates;
	
public:

	/// Blank Constructor
	Cell() {}

	/// insert pointer
	void insert_element(std::shared_ptr<P> _p) { aggregates.push_back(_p); }

	/// Returns the list in the cell
	std::vector<std::shared_ptr<P>>& get_list() { return aggregates; }

	/// Returns number of aggregates in the cell
	int get_n() const { return aggregates.size(); }

	/// Returns the nth element in the cell
	/// \param _idx aggregate index in the cell 
	std::shared_ptr<P>& operator[](const int _idx) { return aggregates[_idx]; };

	/// Erases the nth element in the cell
	void erase_n(int _index) {
		aggregates.erase(aggregates.begin() + _index);
	}

	/// Clears the cell
	void clear_cell() { aggregates.clear(); }

};

#endif
