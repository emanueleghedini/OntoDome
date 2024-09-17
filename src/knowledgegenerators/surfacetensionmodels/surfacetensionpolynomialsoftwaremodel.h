#ifndef SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H
#define SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H

#include "surfacetensionmaterialrelation.h"

class SurfaceTensionPolynomialSoftwareModel : public SoftwareModel {
private:
    /// Core implementation of the model.
    /// It stands fot the numerical representation of the model.
    /// \param T species temperature [K].
    double impl(double T)
    {
        return (s[0] - (s[1]*(T - s[2])));
    }

    std::vector<double> s; ///< Model coefficients array container.
    double* sval; ///< Pointer to the species's property data type.

public:
    SurfaceTensionPolynomialSoftwareModel() : SoftwareModel() {}

    std::string getClassName() const { return "SurfaceTensionPolynomialSoftwareModel"; }

    /// Returns a brief description of the software model by means of its mathematical representation
    /// and the used parameters.
    std::string model_description() {
      return "a - b * ( T - c ) \nwhere a,b and c are the coefficients and T is the temperature.";
    }

    /// Select the model coefficients based on Species Symbol and quantity pointer.
    /// In order to save computational time on several calls of the run method, this
    /// is done on first call only.
    /// Then gathers all the necessary inputs through the relations graph and
    /// runs the actual implementation ("impl").
    void run() {
      if (s.size() == 0) {
        std::string _name = findNearest<SingleComponentComposition>()->name;
        s = get_coeffs(_name);
        sval = &findNearest<SurfaceTension>()->getRelatedObjects<Real>()[0]->value;
      }

      //Get the object's first related temperature
      double T = findNearest<Temperature>()->getRelatedObjects<Real>()[0]->value;

      // Compute the value and push it to the Species' Surface Tension object
      *sval = impl(T);
    }

    /// Return the coefficients needed by the model on given chemical element name.
    /// \param name species name.
    std::vector<double> get_coeffs(std::string _name) const {
      if (_name == "Si") {
        return {0.732, 0.000086, 1685.};
      }
      else if (_name == "Ar") {
        return {0.,0.,0.};
      }
      else if (_name == "H2") {
        return {0.,0.,0.};
      }
      else if (_name == "N2") {
        return {0.,0.,0.};
      }
      else if (_name == "O2") {
        return {0.,0.,0.};
      }
      else if (_name == "Fe") {
        return {1.92,-0.000397,1811.};
      }
      else if (_name == "Ti") {
        return {1.557,0.000156,1941.};
      }
      else if (_name == "Cu") {
        return {1.334,-0.000256,1357.77};
      }
      else if (_name == "Al") {
        return {0.88,-0.000128,933.};
      }
      else if (_name == "He") {
        return {0.,0.,0.};
      }
      else if (_name == "Ag") {
        return {0.894,-0.000191,1234.96};
      }
      else { std::cout << "Species " << _name << " not found!" << std::endl;
          abort(); }
    }

};

#endif // SURFACETENSIONPOLYNOMIALSOFTWAREMODEL_H
