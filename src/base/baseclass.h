#ifndef BASECLASS_H
#define BASECLASS_H

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

boost::uuids::random_generator uuid_gen;

/// Ontology's root class.
/// It provides a UUID for each object when instantiated.

class BaseClass
{
    boost::uuids::uuid uuid; // object's uuid

public:
    BaseClass() { uuid = uuid_gen(); } // every object gets a UUID when instantiated

    /// Returns the object's UUID.
    boost::uuids::uuid getUuid() const { return uuid; }

};

#endif // BASECLASS_H
