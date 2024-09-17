#ifndef SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H
#define SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H

#include "saturationpressurematerialrelation.h"

class SaturationPressurePolynomialSoftwareModel : public SoftwareModel {
private:
    /// Core implementation of the model.
    /// It stands fot the numerical representation of the model.
    ///  \param T species temperature [K].
    double impl(double T)
    {
        return 1.01e5 * pow(10.0,(s[0]-(s[1]/T)));
    }

    std::vector<double> s; ///< Model coefficients array container.
    double* sval; ///< Pointer to the species's property data type.

public:
    SaturationPressurePolynomialSoftwareModel() : SoftwareModel() {}

    std::string getClassName() const { return "SaturationPressurePolynomialSoftwareModel"; }

    /// Returns a brief description of the software model by means of its mathematical representation
    /// and the used parameters.
    std::string model_description() {
      return "10^(a - b/T) \nwhere a and b are the coefficients and T is the temperature.";
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

        sval = findNearest<SaturationPressure>()->value;
      }

      // Get the object's first related temperature
      double T = *findNearest<Temperature>()->get_data();

      // Compute the value and push it to the Species' Saturation Pressure object
      *sval = impl(T);
    }

    /// Return the coefficients needed by the model on given chemical element name.
    /// \param name species name.
    std::vector<double> get_coeffs(std::string _name) const {
      if (_name == "Si") {
        return {7.5341, 23399.};
      }
      else if (_name == "Ar") {
        return {0., 0.};
      }
      else if (_name == "H2") {
        return {0., 0.};
      }
      else if (_name == "N2") {
        return {0., 0.};
      }
      else if (_name == "O2") {
        return {0., 0.};
      }
      else if (_name == "Fe") {
        return {12.535, 21400.};
      }
      else if (_name == "Ti") {
        return {15.5, 55200.};
      }
      else if (_name == "Cu") {
        return {10.845, 16580.};
      }
      else if (_name == "Al") {
        return {5.911, 16211.};
      }
      else if (_name == "He") {
        return {0.,0.};
      }
      else if (_name == "Ag") {
        return {10.785, 13827.};
      }
      else { std::cout << "Species " << _name << " not found!" << std::endl;
          abort(); }
    }
};

#endif // SATURATIONPRESSUREPOLYNOMIALSOFTWAREMODEL_H
