#ifndef THING_H
#define THING_H

#include <iostream>
#include <vector>
#include <map>

#include "baseclass.h"
#include "relation.h"
#include "datatypes.h"
#include "algorithm"

class Thing : public BaseClass {

    // dynamic relations
    std::vector<Relation*> relations;

    // Find related objects through the entity relations hierarchy with class as given
    // uses a list of the already scanned objects' uuids
    template<class T> T* looker(Thing* start, std::vector<boost::uuids::uuid>* scanned);

public:
    virtual std::string getClassName() const { return "Thing"; }

    // populate with relations
    void addRelation(Relation* t) { relations.push_back(t); }

    // removes a relation
    void removeRelation(Relation* r) {
      int i = 0;
      for (auto i : relations) {
          if (r == i) {
            break;
          }
        i += 1;
      }
      relations.erase(relations.begin()+i);
    }

    // Create a relation to this entity
    template<class T, class T0>
    void createRelationTo(T0* o1) {

        T* r = new T(this,o1);

        this->addRelation(r);
        o1->addRelation(r);
    }

    // Creates multiple relations to this entity at once
    template<class T, class T0>
    void createRelationsTo(std::vector<T0*> o1) {

      for (auto i : o1) {
        T* r = new T(this,i);

        this->addRelation(r);
        i->addRelation(r);
      }
    }

    // find a relation of a specific type
    // TODO: optimize search using multimap<std::string,Relation*>,
    //       getClassName() function for the key and a static_cast<T*>

    // Gets all the entities with given relation to this entity
    template<class T> std::vector<T*> getRelations();

    // Gets the last entity with given relation to this entity
    template<class T> T* getLastRelation();

    // Find first related object through the entity relations hierarchy with class as given
    template<class T> T* find();

    // Find all related objects through the entity relations hierarchy with class as given
    template<class T> std::vector<T*> findAll();

    // virtual run method for objects which requires it
    virtual void run() { std::cout << "what to run? eh?" << std::endl;};

    // Gets the all the entities with the specified class related to given entity
    template<class T>
    std::vector<T*> getRelatedObjects();
};

class EMMO : public Thing {
public:
    std::string getClassName() const { return "EMMO"; }
};

class Item : public EMMO {
public:
    std::string getClassName() const { return "Item"; }
};

class Collection : public EMMO {
public:
    std::string getClassName() const { return "Collection"; }
};

class Quantum : public Item {
public:
    std::string getClassName() const { return "Quantum"; }
};

class Void : public Item {
public:
    std::string getClassName() const { return "Void"; }
};

class Physical : public Item {
public:
    std::string getClassName() const { return "Physical"; }
};

class Perspective :virtual public Physical {
public:
    std::string getClassName() const { return "Perspective"; }
};

class Reductionistic : public Perspective {
public:
    std::string getClassName() const { return "Reductionistic"; }
};

class Existent : public Reductionistic {
public:
    std::string getClassName() const { return "Existent"; }
};

class State : public Existent {
public:
    std::string getClassName() const { return "State"; }
};

class ChemicalComposition : public State {
public:
    std::string getClassName() const { return "ChemicalComposition"; }
};

class Equation : public State {
public:
    std::string getClassName() const { return "Equation"; }
};

class MaterialRelation : public Equation {
public:
    std::string getClassName() const { return "MaterialRelation"; }
};

class Symbolic : public Perspective {
public:
    std::string getClassName() const { return "Symbolic"; }
};

class SymbolicConstruct : public Symbolic {
public:
    std::string getClassName() const { return "SymbolicConstruct"; }
};

class String : public SymbolicConstruct, public DataType<std::string>
{
public:
    std::string getClassName() const { return "String"; }
};

class Arrangement : public State
{
public:
    std::string getClassName() const { return "Arrangement"; }
};

class Array : public Arrangement
{
public:
    std::string getClassName() const { return "Array"; }
};

class Vector : public Array, public DataType<std::vector<double>>
{
public:
    Vector(std::vector<double> s) {data = s;}
    std::string getClassName() const { return "Vector"; }
};

class Language : public Symbolic {
public:
    std::string getClassName() const { return "Language"; }
};

class Chemical : public Language {
public:
    std::string getClassName() const { return "Chemical"; }
};

class ChemicalSpecies : public Chemical {
public:
    std::string getClassName() const { return "ChemicalSpecies"; }
};

class ChemicalElement : public ChemicalSpecies {
public:
    std::string getClassName() const { return "ChemicalElement"; }
};

class Metrological : public Language {
public:
    std::string getClassName() const { return "Metrological"; }
};

class IUPAC : public String
{
public:
    IUPAC(std::string _iupac) { data = _iupac; }
    std::string getClassName() const { return "IUPAC Name"; }
};

class LatexExpression : public String
{
public:
    LatexExpression(std::string _expr) { data = _expr; }
    std::string getClassName() const { return "LatexExpression"; }
};

class Unit : public String
{
public:
    Unit(std::string _unit) {data = _unit;}
    std::string getClassName() const { return "Unit"; }
};

class Mathematical : public Language {
public:
    std::string getClassName() const { return "Mathematical"; }
};

class Numerical : public Mathematical {
public:
    std::string getClassName() const { return "Numerical"; }
};

class Number : public Numerical {
public:
    std::string getClassName() const { return "Number"; }
};

class Real : public Number, public DataType<double>
{
public:
    Real(double s) {data = s;}
    std::string getClassName() const { return "Real"; }
};

class KnowledgeGenerator : public Perspective {
public:
    std::string getClassName() const { return "KnowledgeGenerator"; }
};

class Holistic : public Perspective {
public:
    std::string getClassName() const { return "Holistic"; }
};

class Participant : public Holistic {
public:
    std::string getClassName() const { return "Participant"; }
};

class Semiotic : public Participant {
public:
    std::string getClassName() const { return "Semiotic"; }
};

class Sign : public Semiotic {
public:
    std::string getClassName() const { return "Sign"; }
};

class Icon : public Sign {
public:
    std::string getClassName() const { return "Icon"; }
};

class Model : public Icon {
public:
    std::string getClassName() const { return "Model"; }
};

class SoftwareModel : public Model, public KnowledgeGenerator {
public:
    std::string getClassName() const { return "SoftwareModel"; }

    // Calls the model's/material relation's software implementation method.
    // If multiple methods are defined for the same SoftwareModel object, the first added is selected;
    void run() {
      find<SoftwareModel>()->run();
    }
};

class MathematicalModel : public Model {
public:
    std::string getClassName() const { return "MathematicalModel"; }
};

class PhysicsBasedModel : public MathematicalModel {
public:
    std::string getClassName() const { return "PhysicsBasedModel"; }
};

class ContinuumModel : public PhysicsBasedModel {
public:
    std::string getClassName() const { return "ContinuumModel"; }
};

class MesoscopicModel : public PhysicsBasedModel {
public:
    std::string getClassName() const { return "MesoscopicModel"; }
};

class Physicalistic : public Perspective {
public:
    std::string getClassName() const { return "Physicalistic"; }
};

class Matter : public Physicalistic {
public:
    std::string getClassName() const { return "Matter"; }
};

class Continuum : public Matter {
public:
    std::string getClassName() const { return "Continuum"; }
};

class Fluid : public Continuum {
public:
    std::string getClassName() const { return "Fluid"; }
};

class Gas : public Fluid {
public:
    std::string getClassName() const { return "Gas"; }
};

class GasMixture : public Gas {
public:
    std::string getClassName() const { return "GasMixture"; }
};

class ChemicalEntity : public Matter {
public:
    std::string getClassName() const { return "ChemicalEntity"; }
};

class MolecularEntity : public ChemicalEntity {
public:
    std::string getClassName() const { return "MolecularEntity"; }
};

class Atom : public MolecularEntity {
public:
    std::string getClassName() const { return "Atom"; }
};

class PolyatomicEntity : public MolecularEntity {
public:
    std::string getClassName() const { return "PolyatomicEntity"; }
};

class Molecule : public PolyatomicEntity {
public:
    std::string getClassName() const { return "Molecule"; }
};

class HeteronuclearMolecule : public Molecule {
public:
    std::string getClassName() const { return "HeteronuclearMolecule"; }
};

class HomonuclearMolecule : public Molecule {
public:
    std::string getClassName() const { return "HomonuclearMolecule"; }
};

class Quantity : public State, public Metrological
{
public:
    Quantity(Real* _s, Unit* _u)
    {
        State::createRelationsTo<hasPart,Thing>({_u,_s});
    }

    std::string getClassName() const { return "Quantity"; }
};

class PhysicalQuantity : public Quantity
{
public:
    PhysicalQuantity(Real* _s, Unit* _u) : Quantity(_s,_u)
    {
        State::createRelationsTo<hasPart,Thing>({_u,_s});
    }

    std::string getClassName() const { return "PhysicalQuantity"; }
};

class BaseQuantity : public PhysicalQuantity
{
public:
    BaseQuantity(Real* _s, Unit* _u) : PhysicalQuantity(_s,_u)
    {
        State::createRelationsTo<hasPart,Thing>({_u,_s});
    }

    std::string getClassName() const { return "BaseQuantity"; }
};

class Time : public BaseQuantity {
public:
    Time(Real* _s, Unit* _u) : BaseQuantity(_s,_u) {}

    std::string getClassName() const { return "Time"; }
};

class DerivedQuantity : public PhysicalQuantity
{
public:
    DerivedQuantity(Real* _s, Unit* _u) : PhysicalQuantity(_s,_u)
    {
        State::createRelationsTo<hasPart,Thing>({_u,_s});
    }

    std::string getClassName() const { return "DerivedQuantity"; }
};

class Pressure : public DerivedQuantity
{
public:
    Pressure(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "Pressure"; }
};

class Temperature : public BaseQuantity
{
public:
    Temperature(Real* _s, Unit* _u) : BaseQuantity(_s,_u) {}

    std::string getClassName() const { return "Temperature"; }
};

class Mass : public BaseQuantity
{
public:
    Mass(Real* _s, Unit* _u) : BaseQuantity(_s,_u) {}

    std::string getClassName() const { return "Mass"; }
};

class Density : public DerivedQuantity
{
public:
    Density(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "Density"; }
};

class BulkDensityLiquid : public Density
{
public:
    BulkDensityLiquid(Real* _s, Unit* _u) : Density(_s,_u) {}

    std::string getClassName() const { return "BulkDensityLiquid"; }
};

class BulkDensitySolid : public Density
{
public:
    BulkDensitySolid(Real* _s, Unit* _u) : Density(_s,_u) {}

    std::string getClassName() const { return "BulkDensitySolid"; }
};

class ISQDimensionlessQuantity : public DerivedQuantity
{
public:
    ISQDimensionlessQuantity(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "ISQDimensionlessQuantity"; }
};

class RatioQuantity : public ISQDimensionlessQuantity
{
public:
    RatioQuantity(Real* _s, Unit* _u) : ISQDimensionlessQuantity(_s,_u) {}

    std::string getClassName() const { return "RatioQuantity"; }
};

class MolarFraction : public RatioQuantity
{
public:
    MolarFraction(Real* _s, Unit* _u) : RatioQuantity(_s,_u) {}

    std::string getClassName() const { return "MolarFraction"; }
};

class Viscosity : public DerivedQuantity
{
public:
    Viscosity(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "Viscosity"; }
};

class DynamicViscosity : public Viscosity
{
public:
    DynamicViscosity(Real* _s, Unit* _u) : Viscosity(_s,_u) {}

    std::string getClassName() const { return "DynamicViscosity"; }
};

class KinematicViscosity : public Viscosity
{
public:
    KinematicViscosity(Real* _s, Unit* _u) : Viscosity(_s,_u) {}

    std::string getClassName() const { return "KinematicViscosity"; }
};

class SurfaceTension : public DerivedQuantity
{
public:
    SurfaceTension(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "SurfaceTension"; }
};

class SaturationPressure : public DerivedQuantity
{
public:
    SaturationPressure(Real* _s, Unit* _u) : DerivedQuantity(_s,_u) {}

    std::string getClassName() const { return "SaturationPressure"; }
};

class SingleComponentComposition : public ChemicalComposition {
public:
    MolarFraction* mol;
    std::string name;

    SingleComponentComposition(MolarFraction* _mol, std::string _name) { mol = _mol; name = _name; }
    std::string getClassName() const { return "SingleComponentComposition"; }
};

#include "thing.cpp"

#include "species.h"

#endif // THING_H
