#ifndef RELATION_H
#define RELATION_H

#include <string>
#include <vector>

#include "baseclass.h"

/// This class includes all the ontological relations that can be defined
/// between two objects. They follow the hierarchical structure defined under the EMMO project.

class Thing;

class Relation : public BaseClass {
    Thing *o0;
    Thing *o1;

public:
    Relation(Thing* _o0, Thing* _o1)
        : o0(_o0), o1(_o1) {}

    /// Returns the name of the defined relation.
    virtual std::string getRelationName() const { return "Relation"; }

    /// Returns the Domain of a relation (the object for which the relation was defined).
    Thing* getDomain() { return o0; }

    /// Returns the Range of a relation (the object to which the relation points).
    Thing* getRange()  { return o1; }

    virtual ~Relation() = default;
};

class hasOverlap : public Relation {
public:
    hasOverlap(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}
};

class hasPart : public hasOverlap {
public:
    hasPart(Thing* _o0, Thing* _o1) : hasOverlap(_o0,_o1) {}

    virtual std::string getRelationName() const { return "hasPart"; }
};

class hasProperPart : public hasPart {
public:
    hasProperPart(Thing* _o0, Thing* _o1) : hasPart(_o0,_o1) {
        if(_o0==_o1) { abort(); }
    }

    std::string getRelationName() const { return "hasProperPart"; }
};

class hasTemporalPart : public hasProperPart {
public:
    hasTemporalPart(Thing* _o0, Thing* _o1) : hasProperPart(_o0,_o1) {}

    virtual std::string getRelationName() const { return "hasTemporalPart"; }
};

class hasTemporalDirectPart : public hasTemporalPart {
public:
    hasTemporalDirectPart(Thing* _o0, Thing* _o1) : hasTemporalPart(_o0,_o1) {}

    virtual std::string getRelationName() const { return "hasTemporalDirectPart"; }
};

class hasSpatialPart : public hasProperPart {
public:
    hasSpatialPart(Thing* _o0, Thing* _o1) : hasProperPart(_o0,_o1) {}

    virtual std::string getRelationName() const { return "hasSpatialPart"; }
};

class hasConnection : public Relation {
public:
  hasConnection(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}

  std::string getRelationName() const { return "hasConnection"; }
};

class hasSign : public Relation {
public:
  hasSign(Thing* _o0, Thing* _o1) : Relation(_o0,_o1) {}
};

class hasInput : public hasProperPart {
public:
    hasInput(Thing* _o0, Thing* _o1) : hasProperPart(_o0,_o1) {}
};

class hasOutput : public hasProperPart {
public:
    hasOutput(Thing* _o0, Thing* _o1) : hasProperPart(_o0,_o1) {}
};

class hasModel : public hasSign {
public:
    hasModel(Thing* _o0, Thing* _o1) : hasSign(_o0,_o1) {}
};

class hasMathematicalModel : public hasModel {
public:
    hasMathematicalModel(Thing* _o0, Thing* _o1) : hasModel(_o0,_o1) {}
};

class hasSoftwareModel : public hasModel {
public:
    hasSoftwareModel(Thing* _o0, Thing* _o1) : hasModel(_o0,_o1) {}
};

class hasProperty : public hasSign {
public:
  hasProperty(Thing* _o0, Thing* _o1) : hasSign(_o0,_o1) {}

  std::string getRelationName() const { return "hasProperty"; }
};

class hasScalarProperty : public hasProperty {
public:
  hasScalarProperty(Thing* _o0, Thing* _o1) : hasProperty(_o0,_o1) {}

  std::string getRelationName() const { return "hasScalarProperty"; }
};


// Relation based functions
template<class T, class T0, class T1>
Relation* addRelation(T0* o0, T1* o1);

#include "relation.cpp"

#endif // RELATION_H
