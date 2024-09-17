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

#ifndef GASMODEL_H
#define GASMODEL_H

#include <valarray>
#include <map>

#include "gasmodels.h"

/// Class implementing the gas phase.
/// The gas phase is univocally determined knowing pressure, temperature and species molar fractions.
/// Access to species properties is provided by hasPart relations defined in the GasMixture object.
class GasModel  : public GasModels {

protected:
    Pressure* p; ///< Pointer to gas phase pressure [Pa].
    Temperature* T; ///< Pointer to gas phase temperature [K].
    PressureTimeDerivative* dpdt; ///< Pointer to gas phase pressure time derivative [Pa/s].
    TemperatureTimeDerivative* dTdt; ///< Pointer to gas phase temperature time derivative [K/s].
    double gamma; ///< expansion coefficient [1/s].
    std::vector<SingleComponentComposition*> specs; ///< vector containing all the species.

public:
    GasModel() {}

    std::string getClassName() const { return "Gas Continuum Model"; }

protected:
    /// Initialize the model by looking for the required inputs through the relations graph.
    void initialize() {

      // Get the GasMixture conditions and components
      GasMixture* gp = findNearest<GasMixture>();
      specs = gp->getRelatedObjects<SingleComponentComposition>();
      p = gp->findNearest<Pressure>();
      dpdt = gp->findNearest<PressureTimeDerivative>();
      T = gp->findNearest<Temperature>();
      dTdt = gp->findNearest<TemperatureTimeDerivative>();

      // Get the species molar fractions and Populate the hash map for name to index resolution
      if (!specs.empty())
      {
          w.resize(specs.size());
          for (std::size_t i = 0; i < specs.size(); ++i)
          {
            w[i] = specs[i]->mol->get_data();
            hash[specs[i]->name] = i;
          }
      }
      else
      {
        abort();
      }

      // Mark the model as initialized
      init = true;
    }

public:
    /// Run the model for a given timestep and species consumption array.
    /// \param dt temporal step [s].
    /// \param w_cons species consumption array [#/m3/s].
    void timestep(double dt, std::valarray<double> w_cons) {
      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Check if molar fraction vector total sum is smaller than 1 for consistency
      // and re-normalize.
      double w_sum = 0;
      for (std::size_t i = 0; i < w.size(); ++i) { w_sum += *w[i]; }
      if (w_sum >= 1.) {
          // std::cout << "Sum of species molar fractions greater than 1.0. Aborting." << std::endl;
          // abort();
          for (std::size_t i = 0; i < w.size(); ++i) { *w[i] /= w_sum; }
      }

      // Check if consumption vector and available species number matches
      // If not, the code must abort
      if (specs.size() != w_cons.size()) { abort(); }

      // Compute the new species number densities
      double n = get_n(); // number of monomers on m3
      double w_cons_tot = w_cons.sum(); // total number of consumed monomers

      std::valarray<double> ns(w.size()); // number of monomers for each species
      for (std::size_t i = 0; i < w.size(); ++i) { ns[i] = *w[i] * n; }

      gamma = w_cons_tot/n + *dTdt->get_data() / *T->get_data() - *dpdt->get_data() / *p->get_data();

      // simple explicit ODE timestep
      // equation is solved for the number density
      for (std::size_t i = 0; i < w_cons.size(); ++i) {
          ns[i] += (w_cons[i] - ns[i] * gamma) * dt; // Volume contraction
          if (ns[i] < 0) {
            ns[i] = 0;
          }
      }

      // gas phase pressure and temperature update
      T->set_data(*T->value + *dTdt->get_data() * dt);
      p->set_data(*p->value + *dpdt->get_data() * dt);

      // molar fractions update
      n = get_n();
      for (std::size_t i = 0; i < w.size(); ++i) { *w[i] = ns[i] / n; }
    }

    /// Update c composition for reactor network
    void c_update() {
      for (std::size_t i = 0; i < w.size(); ++i) {
        c[i] = *w[i];
      }
    }

    /// Updates gas' state
    void update(double _p, double _T, std::valarray<double> _c) {
      // Initialize the model if not done before
      if (init == false) { initialize(); }

      p->set_data(_p);
      T->set_data(_T);
      c = _c;
    }

    /// Get the expansion coefficient [1/s].
    double get_gamma() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return gamma; }

    /// Get the temperature [K].
    double get_T() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return *T->value; }

    /// Get the pressure [pa].
    double get_p() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return *p->value; }

    /// Get species molar fractions [%]
    std::valarray<double> get_c() { return c; }

    /// Get gas phase molecules average viscosity [Pa s].
    double get_average_viscosity() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      double visc = 0;

      for(size_t i=0; i<w.size(); ++i)
       visc += *specs[i]->getRelatedObjects<Viscosity>()[0]->get_data() * *w[i];

      return visc;
    }

    /// Get gas phase number density [#/m3].
    double get_n() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return get_p()/(K_BOL* get_T());
    }

    /// Get the superaturation ratio [#].
    /// \param spec selected species.
    template <class TT> double get_S(TT* spec) {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      double m_frac = *w[hash.at(spec->name)];
      double ns = m_frac*get_n();
      double n_sat = get_n_sat<TT>(spec);

      return ns/n_sat;
    }

    /// Get the saturation density [#/m3].
    /// \param spec selected species.
    template <class TT> double get_n_sat(TT* spec) {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      // Update specie's saturation pressure
      spec->template getRelatedObjects<SaturationPressure>()[0]->template getRelatedObjects<SaturationPressureMaterialRelation>()[0]->run();

      return *spec->template getRelatedObjects<SaturationPressure>()[0]->template value /(K_BOL*get_T());
    }

    /// Get the gas phase mass density [kg/m3].
    double get_density() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return get_average_molecular_mass() * get_p()/(K_BOL*get_T()); }

    /// Get the gas phase molecules average mass [kg].
    double get_average_molecular_mass() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      double m = 0.;

      for(size_t i=0; i<w.size(); ++i)
          m += *specs[i]->getRelatedObjects<Mass>()[0]->get_data() * *w[i];

      return m;
    }

    /// Get the gas phase mean free path [m].
    double get_mfp() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      return (get_average_viscosity()/get_p()) * sqrt(M_PI*K_BOL*get_T()/(2.*get_average_molecular_mass())); }

    /// Get the gas flux used for Langevin dynamics [kg/m2 s].
    double get_gas_flux() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      double flux = 0.;

      for(size_t i=0; i<w.size(); ++i) {

          double m_gas = *specs[i]->getRelatedObjects<Mass>()[0]->get_data();

          // n_s * m_s * sqrt(3*K_BOL*T/m_s);
          flux += *w[i]*get_p()/(K_BOL*get_T()) * m_gas * sqrt(3*K_BOL*get_T()/m_gas);
      }
      return flux;
    }

    /// Prints the gas phase most relevant properties.
    void print() {

      // Initialize the model if not done before
      if (init == false) { initialize(); }

      std::cout << "Time: " << *findNearest<Time>()->get_data() << std::endl;

      for(auto& sp : specs)
          std::cout << sp->name << '\t';
      std::cout << std::endl;

      for(auto& cs : w)
          std::cout << *cs << '\t';
      std::cout << std::endl << std::endl;

      std::cout << "T[K]\t" << "p[Pa]\t" << "n[#/m3]\t  " << "gamma[1/s]\t" << std::endl;
      std::cout << get_T() << "\t" << get_p() << "\t" << get_n() << "\t" << "  " << get_gamma() << "\t" << std::endl << std::endl;
    };

};

#endif // GASMODEL_H
