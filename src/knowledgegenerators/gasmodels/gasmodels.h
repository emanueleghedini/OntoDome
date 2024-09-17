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

#ifndef GASMODELS_H
#define GASMODELS_H

#include "../../ontodome.h"

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by hasPart relations defined in the GasMixture object.
class GasModels  : public SoftwareModel {

protected:

    std::valarray<double*> w; ///< Species molar fractions.
    std::valarray<double> c; ///< Species molar fractions for reactor network.
    std::map<std::string,std::size_t> hash; ///< Hash map for name-to-index resolution.
    bool init = false; ///< Boolean value which stores whether the model is initialized or not.

public:

    GasModels() {}

    std::string getClassName() const { return "Gas Continuum Models Container"; }

    /// Initialize the model by looking for the required inputs through the relations graph.
    virtual void initialize() = 0;

    /// Run the model for a given timestep and species consumption array.
    /// \param dt temporal step [s].
    /// \param w_cons species consumption array [#/m3/s].
    virtual void timestep(double dt, std::valarray<double> w_cons) = 0;

    /// Get the expansion coefficient [1/s].
    virtual double get_gamma() = 0;

    /// Get the temperature [K].
    virtual double get_T() = 0;

    /// Get the pressure [pa].
    virtual double get_p() = 0;

    /// Get the nanonetwork molar fractions [#].
    virtual std::valarray<double> get_c() = 0;

    /// Get gas phase molecules average viscosity [Pa s].
    virtual double get_average_viscosity() = 0;

    /// Get gas phase number density [#/m3].
    virtual double get_n() = 0;

    /// Get the superaturation ratio [#].
    /// \param spec selected species.
    template <class TT> double get_S(TT* spec) {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      double m_frac = *w[hash.at(spec->name)];
      double ns = m_frac*get_n();
      double n_sat = get_n_sat<TT>(spec);

      return ns/n_sat;
    };

    /// Get the saturation density [#/m3].
    /// \param spec selected species.
    template <class TT> double get_n_sat(TT* spec) {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Update specie's saturation pressure
      spec->template getRelatedObjects<SaturationPressure>()[0]->template getRelatedObjects<SaturationPressureMaterialRelation>()[0]->run();

      return *spec->template getRelatedObjects<SaturationPressure>()[0]->value /(K_BOL*get_T());
    };

    /// Get the gas phase mass density [kg/m3].
    virtual double get_density() = 0;

    /// Get the gas phase molecules average mass [kg].
    virtual double get_average_molecular_mass() = 0;

    /// Get the gas phase mean free path [m].
    virtual double get_mfp() = 0;

    /// Get the gas flux used for Langevin dynamics [kg/m2 s].
    virtual double get_gas_flux() = 0;

    /// Prints the gas phase most relevant properties.
    virtual void print() = 0;

    /// Updates gas phase pressure, temperature and molar fractions
    /// \param p pressure [Pa]
    /// \param T temperature [K]
    /// \param c valarray with molar fraction of each species ordered as _species
    virtual void update(double p, double T, std::valarray<double> c) = 0;

    /// Updates gas phase molar fractions for reactor network
    virtual void c_update() = 0;

};

#endif // GASMODEL_H
