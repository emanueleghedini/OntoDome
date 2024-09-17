#ifndef DATATYPES_H
#define DATATYPES_H

/// Generic data type base template class.
/// It provides the data entry to be filled with any data
/// such as double, float, string etc.

template<class T>
struct DataType
{
    T value;

    DataType() {}
    DataType(T _value) : value(_value) {}
};

#endif // DATATYPES_H
