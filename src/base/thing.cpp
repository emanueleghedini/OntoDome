#include "thing.h"

// internal tools
template <typename T>
bool contains(std::vector<T> vec, const T & elem)
{
    return any_of(vec.begin(), vec.end(), [&](const auto & x){
        return x == elem;
    });
}

// sources
std::vector<Thing *> Thing::ListRelations() {
  std::vector<Thing *> res;
    for (auto x : this->relations) {
        if (x->getRange()->getUuid() != this->getUuid()) {
            res.push_back(x->getRange());
        } else {
            res.push_back(x->getDomain());
        }
    }
   return res;
}

void Thing::addRelation(Relation* t) {
  if (std::find(relations.begin(), relations.end(), t) == relations.end()) {
  relations.push_back(t);
  } else { abort(); }
}

void Thing::removeRelation(Relation* r) {
  int i = 0;
  // look for the index associated to the relation
  for (auto i : relations) {
    if (r == i) { break; }
    i += 1;
  }
  if (i != 0) { relations.erase(relations.begin()+i); }
  else { abort(); }
}

template<class T, class T0>
void Thing::createRelationTo(T0* o1) {

    T* r = new T(this,o1);

    this->addRelation(r);
    o1->addRelation(r);
}

template<class T, class T0>
void Thing::createRelationsTo(std::vector<T0*> o1) {

  for (auto i : o1) {
    T* r = new T(this,i);

    this->addRelation(r);
    i->addRelation(r);
  }
}

template<class T>
std::vector<T*> Thing::getRelations() {
   std::vector<T*> rels;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(i)) {
        rels.push_back(r0);
      }
   }

   return rels;
}

/// Avoids looping on same object by returning the other object in a given relation.
Thing* avoider (Thing* o0, Relation* o1) {
  Thing* tested;
  if (o1->getRange()->getUuid() != o0->getUuid()) {
      tested = o1->getRange();
  } else {
      tested = o1->getDomain();
  }
  return tested;
}

/// Adds the UUID of an object only if not already present in the list.
void adder(std::vector<boost::uuids::uuid>* scans, boost::uuids::uuid id) {
  if (std::find(scans->begin(), scans->end(), id) == scans->end()) {
    scans->push_back(id);
  }
}

template<class T>
T* Thing::looker(Thing* start, std::vector<boost::uuids::uuid>* scanned) {
  T* result = nullptr;
  adder(scanned,start->getUuid());
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    if (std::find(scanned->begin(), scanned->end(), tested->getUuid()) == scanned->end()) {
      adder(scanned,tested->getUuid()); // add to scanned list
      if (dynamic_cast<T*>(tested)) {
        result = dynamic_cast<T*>(tested);
        return result;
        break;
      }
    }
  }

  if (result == nullptr) {
    for (auto i : start->relations) {
      auto tested = avoider(start,i);
      adder(scanned,tested->getUuid()); // add to scanned list

      // look in tested object's relations, if there is an unscanned object, mark tested to be looked.
      bool to_scan = false;
      for (auto ii : tested->relations) {
          if (std::find(scanned->begin(), scanned->end(), avoider(tested,ii)->getUuid()) == scanned->end()) {
              to_scan = true;
              break;
          }
      }

      if (to_scan) {
        result = looker<T>(tested, scanned);
        if (result != nullptr) { break; }
      }
    }
  }
  return result;
}

template<class T>
T* Thing::findNearest() {
  auto* start = this; // the starting object will always be the objects from which the method is called
  T* result = nullptr;
  std::vector<boost::uuids::uuid> scanned; // list of scanned objects by means of UUID
  adder(&scanned,start->getUuid()); // add start to scanned items to avoid looping

  // scan for object of type T in start object relations
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); //add to scanned list
    if (dynamic_cast<T*>(tested)) {
      result = dynamic_cast<T*>(tested);
      break; // stop searching, we got it
    }
  }

  // if not found above, try searching in related objects
  if (result == nullptr) {
    for (auto i : start->relations) {
      auto tested = avoider(start,i);
      adder(&scanned,tested->getUuid()); // add to scanned list
      result = looker<T>(tested, &scanned);
      if (result != nullptr) { break; } // stop searching, we got it
    }
  }

  if (result == nullptr) {
      std::cout << "Not found." << std::endl; // for test purposes
      abort(); }
  else { return result; }
}

template<class T>
std::vector<T*> Thing::findAll() {
  auto* start = this;  // the starting object will always be the objects from which the method is called
  std::vector<T*> result;
  std::vector<boost::uuids::uuid> scanned; // list of scanned objects by means of UUID
  adder(&scanned,start->getUuid()); // add start to scanned items to avoid looping

  // scan for object of type T in start object relations
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); //add to scanned list
    if (dynamic_cast<T*>(tested)) {
      result.push_back(dynamic_cast<T*>(tested));
    }
  }

  // if not found above, try searching in related objects
  for (auto i : start->relations) {
    auto tested = avoider(start,i);
    adder(&scanned,tested->getUuid()); // add to scanned list
    auto res = looker<T>(tested, &scanned);
    if (res != nullptr) {
      result.push_back(looker<T>(tested, &scanned));
    }
  }

  if (result.size() == 0) {
      std::cout << "Nothing found." << std::endl; // for test purposes
      abort(); }
  else { return result; }
}

template<class T>
std::vector<T*> Thing::getRelatedObjects() {
   std::vector<T*> objs;

   for(auto i: relations) {
      if(T* r0 = dynamic_cast<T*>(avoider(this,i))) {
        if (!contains(objs, r0)) {
          objs.push_back(r0);
        }
      }
   }

   return objs;
}
