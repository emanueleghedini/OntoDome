#ifndef SPECIES_H
#define SPECIES_H

/// This class contains all the chemical elements that will be used in ontodome.
/// Each element has a Symbol which corresponds to its standard name (i.e. Periodic Table of Elements)
/// and a method ("get_symbol") which returns it.

#include "thing.h"

class HeliumSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string HeliumSymbol::name = "Helium";

class Helium : public Atom {
public:
    Helium () : Atom () {
      HeliumSymbol* he = 0;
      this->createRelationTo<hasProperty>(he);
    }
    std::string getClassName() const { return "Helium"; }
};

//--------------------------------------------

class ArgonSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string ArgonSymbol::name = "Argon";

class Argon : public Atom {
public:
    Argon () : Atom () {
      ArgonSymbol* ar = 0;
      this->createRelationTo<hasProperty>(ar);
    }
    std::string getClassName() const { return "Argon"; }
};

//--------------------------------------------

class SiliconSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string SiliconSymbol::name = "Silicon";

class Silicon : public Atom {
public:
    Silicon () : Atom () {
      SiliconSymbol* si = 0;
      this->createRelationTo<hasProperty>(si);
    }
    std::string getClassName() const { return "Silicon"; }
};

//--------------------------------------------

class IronSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string IronSymbol::name = "Iron";

class Iron : public Atom {
public:
    Iron () : Atom () {
      IronSymbol* fe = 0;
      this->createRelationTo<hasProperty>(fe);
    }
    std::string getClassName() const { return "Iron"; }
};


//--------------------------------------------

class AluminumSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string AluminumSymbol::name = "Aluminum";

class Aluminum : public Atom {
public:
    Aluminum () : Atom () {
      AluminumSymbol* al = 0;
      this->createRelationTo<hasProperty>(al);
    }
    std::string getClassName() const { return "Aluminum"; }
};

//--------------------------------------------

class SilverSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string SilverSymbol::name = "Silver";

class Silver : public Atom {
public:
    Silver () : Atom () {
      SilverSymbol* ag = 0;
      this->createRelationTo<hasProperty>(ag);
    }
    std::string getClassName() const { return "Silver"; }
};

//--------------------------------------------

class TitaniumSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string TitaniumSymbol::name = "Titanium";

class Titanium : public Atom {
public:
    Titanium () : Atom () {
      TitaniumSymbol* ti = 0;
      this->createRelationTo<hasProperty>(ti);
    }
    std::string getClassName() const { return "Titanium"; }
};

//--------------------------------------------

class CopperSymbol : public ChemicalElement {
  static std::string name;

public:
  static std::string get_symbol() { return name; }
};

std::string CopperSymbol::name = "Copper";

class Copper : public Atom {
public:
    Copper () : Atom () {
      CopperSymbol* cu = 0;
      this->createRelationTo<hasProperty>(cu);
    }
    std::string getClassName() const { return "Copper"; }
};

#endif // SPECIES_H
