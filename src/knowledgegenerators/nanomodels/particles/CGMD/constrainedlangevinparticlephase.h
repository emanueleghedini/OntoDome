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

#ifndef CONSTRAINEDLANGEVINPARTICLEPHASE_H
#define CONSTRAINEDLANGEVINPARTICLEPHASE_H

#include "../particlephase/dynamicparticlephase.h"
#include "../PBM/pbmfractalparticlephase.h"

#include <time.h>

//#define OMP
//#include <omp.h>

#include "../utilities/ndm_random.h"

// External Libs
#include "../utilities/tinyxml2/tinyxml2.h"
#include "../../../../tools/utilities.h"

template<typename A>
class ConstrainedLangevinParticlePhase : public DynamicParticlePhase<A> {

  /// Function for calculating the value of a 3d index a grid in a periodic domain
  ///	\param	int c_idx: cell actual component position
  ///	\param int SIDE: Number of cells in the dimension
  ///	\param double& shift: shift to impose in the point coordinates at the correspondig component
  ///	\param double domain_side: lenght of the control volume side
  int periodic_position(int c_idx, int SIDE, double& shift, double domain_side);

  /// Common variables declaration.
  GasModels* gp; ///< Pointer to currently used GasModel
  NucleationTheory* nt; ///< Pointer to currently used Nucleation theory
  SingleComponentComposition* sp; ///< Pointer to the species for which the PBM method instance is defined
  double* T; ///< Pointer to gas temperature

  bool init = false; ///< Boolean value which stores whether the model is initialized or not.

public:

  /// Constructor
  ConstrainedLangevinParticlePhase(double _volume):
      DynamicParticlePhase<A>(_volume){}

  /// Constructor from PBM Particle phase
  ///	\param PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp PBM particle phase
  ConstrainedLangevinParticlePhase(PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp, double _volume, double T) :
    DynamicParticlePhase<A>(_volume)
  {

    double PBM_Volume = _fractal_pp.volume;

    int N = _fractal_pp.get_aggregates_number();

    //if volumes are the same

    if (PBM_Volume == _volume) {
      // Spawn aggregates from PBM aggregates
      for (int i = 0; i < N; i++) {
        std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(i);
        this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());

      }
    }
    else if (_volume < PBM_Volume) {

      // Choose ramdom aggregates from the PBM Particle Phase
      double pbm_p_density = N / PBM_Volume;
      int n_particles = (N*_volume) / PBM_Volume;
      std::cout << "PBM Particles:" << N << std::endl;
      std::cout << "Langevin Particles" << n_particles << std::endl;

      srand(time(NULL));
      int count = 0;
      while (count < n_particles) {
        int rand_index = rand() % N;
        std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(rand_index);
        this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());
        count++;
      }

    }
    else { // if the volume is bigger than the PBM one
      // TO BE IMPLEMENTED
    }

    // Initialize velocities
    this->initialize_velocities(T);
  }

  /// Constructor for Particle Phase starting from a list of aggregates (hotstart/recovery)
  ConstrainedLangevinParticlePhase(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& _aggregates, double _volume, double _T):
    DynamicParticlePhase<A>(_volume) {

    // Assign the list
    this->aggregates = _aggregates;

    // Initialize velocities
    this->initialize_velocities(_T);

    // Displace aggregates in the grid
    this->grid.displace_aggregates(this->aggregates);
  }

  /// Saves VTK data
  ///	\param int iteration: iteration in which the snapshot is saved
  ///	\param s_idx: int streamline index
  void save_vtk(int iteration, std::string root_path);

  /// Time step
  ///	\param double dt: timestep [sec]
  ///	\param const GasPhase& gp: GasPhase object [by reference]
  ///	\param const NucleationTheory& nt: Nucleation Theory object [by reference]
  ///	\param double& T: Temperature[K] [by reference]
  double timestep(double dt);

  /// Returns cardinality of the active particles in the simulation
  int get_aggregates_cardinality();

  /// Updates graph constraints values
  void update_constraints();

  /// biggest spherical enclosure diameter
  double get_biggest_spherical_enclosure(int& _id) const;


protected:

  /// Initialize the method - attempt to save computational time
  void initialize() {

    // Get all the required inputs

    nt = this->template findNearest<NucleationTheory>();

    gp =  this->template findNearest<GasModels>();

    sp =  this->template findNearest<SingleComponentComposition>();

    T = gp->template findNearest<Temperature>()->get_data();

    init = true;
  }

  /// Motion
  ///	\param dt timestep [s]
  ///	\param gp Gas Phase object
  ///	\param T Temperature [K]
  void motion(double dt, double T);

  ///Coagulation
  void coagulation();

  /// Creates a Particle in the simulation
  void add_particle(double n);

  /// Checks if some particles has coalesced
  void rearrange_aggregate();

  /// Keeps constant the number of aggregates in the simulation
  void aggregates_number_balance();

  /// Initialize Velocity
  ///	\param double T Temperature [K]
  void initialize_velocities(double T);

  /// Evaluate nucleation events for species s and return the molecules consumption
  /// [#/m3 s]
  /// \param dt timestep [s]
  /// \param J particle nucleation rate [#/m3 s]
  /// \param s condensing species
  double nucleation(double dt, double J, double j);

  /// Expands the volume and shifts the aggregates
  /// Correct volume due to gas expansion
  /// \param dt timestep [s]
  /// \param gp gas phase surrounding particles
  void volume_resize(double dt);


};

template<typename A>
int ConstrainedLangevinParticlePhase<A>::periodic_position(int c_idx, int SIDE, double& shift, double domain_side) {
  if (c_idx < 0) { shift = +domain_side; return SIDE - 1; }
  else if (c_idx > (SIDE - 1)) { shift = -domain_side; return 0; }
  else {
    return c_idx;
  }

}

std::vector<std::tuple<int, int, int>> getConnections(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

// vect[0] -> couples
// vect[1] -> types
// vect[2] -> offsets
std::vector<std::vector<int>> GetAbsoluteIndices(std::vector<std::tuple<int, int, int>> graph_data, std::vector<int> p_list);

template<class dataType>std::string V_To_S(std::vector<dataType> vect);

template<typename A>
double ConstrainedLangevinParticlePhase<A>::timestep(double dt) {

  // Initialize the model if not done before
  if (init == false) { initialize(); }

  double g_si = 0.0; // Consumption

  // get Temperature from the Gas Phase
  double Temp = *T;

  // nucleation parameters
  double J_si = nt->nucleation_rate();
  double j_si = nt->stable_cluster_diameter();

  // Condensation parameters
  double F_si = nt->condensation_rate();

  // Consumption of species from the gas phase and nucleation and condensation in the particle phase
  g_si += nucleation(dt, J_si, j_si);

      // Condensation
  g_si += this->condensation(dt, F_si);

  // Move particles
  motion(dt, Temp);

  // Coagulation
  coagulation();

  // Sintering
  this->sintering(dt, Temp);

  // Check Coalescence and applies shake after sintering
  rearrange_aggregate();

  // Keeps constant the number of aggregates
  this->aggregates_number_balance();

  // Update volume from gas phase data
  this->volume_expansion(dt);

  return g_si;
}



template<typename A>
void ConstrainedLangevinParticlePhase<A>::volume_resize(double dt) {

  // Initialize the model if not done before
  if (init == false) { initialize(); }

    double scale_factor = exp(gp->get_gamma()*dt);
    this->volume = this->volume* scale_factor;

        double shift = pow(scale_factor, 1.0 / 3.0);

        std::valarray<double> x_shift = { shift,shift,shift };

    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
                (*agg)->scale_coordinates(x_shift);
        }

}


template<typename A>
void ConstrainedLangevinParticlePhase<A>::rearrange_aggregate() {

  for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
    // Check if the aggregates has coalesced
    if ((*agg)->check_coalescence()) {
      // if yes, update graph structure
      (*agg)->update_graph();
      // Update actual cardinality
      (*agg)->update_particles_cardinality();
    }
    else {
      // Update the constraints with the new distances
      (*agg)->update_constraints();
    }
  }

  std::list<int> erased_agg;

#ifdef OMP

  bool SHAKE_ERROR;

#pragma omp parallel
#pragma omp single nowait
  {
    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
#pragma omp task firstprivate(agg)
      SHAKE_ERROR = (*agg)->SHAKE();
      // Apply shake to the sintered aggregate
#pragma omp critical
      if (!(SHAKE_ERROR)) {
        std::cout << "SHAKE DIVERGENCE IN REARRANGING AGGREGATES AFTER SINTERING" << std::endl;
        erased_agg.push_back((*agg)->get_id());
        this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
      }
    }
#pragma omp taskwait
  }
#else
  for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
    // Apply shake to the sintered aggregate
    if (!((*agg)->SHAKE())) {
      std::cout << "SHAKE DIVERGENCE IN REARRANGING AGGREGATES AFTER SINTERING" << std::endl;
      erased_agg.push_back((*agg)->get_id());
      this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
    }
  }
#endif


    if (erased_agg.size() > 0) {
        // Erase from aggregate
        for(auto e = erased_agg.begin(); e != erased_agg.end(); e++){
            for(auto a = this->aggregates.begin(); a !=this->aggregates.end(); a++){
                if((*a)->get_id() == (*e)){
                    this->aggregates.erase(a);
                    break;
                }
            }
            // Erase from grid
            for(int c = 0; c < this->grid.get_n_cells(); c++){
                for(int a = 0; a < this->grid[c].get_n(); a++){
                    if(this->grid[c][a]->get_id() == (*e)){
                        this->grid[c].erase_n(a);
                        break;
                    }
                }
            }
        }
    }

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::add_particle(double n) {

    // Initialize the model if not done before
    if (this->init == false) { this->initialize(); }

    // Control volume side
    double side = pow(this->volume, 1.0 / 3.0);

    // generate position of the particle
    std::valarray<double> x = { ndm::uniform_double_distr(ndm::rand_gen),
        ndm::uniform_double_distr(ndm::rand_gen),
        ndm::uniform_double_distr(ndm::rand_gen) };

    // scale coordinates on the control volume
    x = side*(x - 0.5);

    // create particle
    std::shared_ptr<DynamicParticle> data(new DynamicParticle(n, this->sp, x));
    std::shared_ptr<RATTLEAggregate<DynamicParticle>> agg(new RATTLEAggregate<DynamicParticle>(data));


    // Add particle to the ensamble.
    this->aggregates.push_back(agg);

    return;

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::initialize_velocities(double T) {

    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
                (*agg)->v_initialization(T);
        }

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::update_constraints() {

    for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
                (*agg)->update_constraints();
        }

}

/// biggest spherical enclosure diameter
template<typename A>
double ConstrainedLangevinParticlePhase<A>::get_biggest_spherical_enclosure(int& _id) const {

  double d;
  double d_max = DBL_MIN;
    for (auto& a : this->aggregates) {
                d = a->get_enclosing_sphere_diameter();
                if (d > d_max) {
                        d_max = d;
                        _id = a->get_id();
                }
        }

        return d_max;
}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::save_vtk(int iter, std::string root_path) {

    using namespace tinyxml2;

    // All points (nanoparticles) coords
    std::vector<double> coords = getParticleCoords(this->aggregates);

    // All points (nanoparticles) mass
    std::vector<double> masses = getParticleMass(this->aggregates);

    // All points (nanoparticles) diameter
    std::vector<double> diameter = getParticleDiameter(this->aggregates);

  // All aggregates diameter
  std::vector<double> agg_diameter = getAggregatesDiameter(this->aggregates);

    // All points (nanoparticles) affiliation
    std::vector<int> affil = getParticleAggregate(this->aggregates);

  // All points IDs
  std::vector<int> p_IDs = getParticleID(this->aggregates);

  // Get the connectivity data
  std::vector <std::tuple<int,int,int>> graph_data = getConnections(this->aggregates);

  std::vector<std::vector<int>> abs_idx = GetAbsoluteIndices(graph_data, p_IDs);

    // Number of elementary nanoparticle in the system
    int nPart = masses.size();

  // Total number of edges in the simulation
  int nEdges = getTotalEdges(this->aggregates);

    XMLDocument xmlDoc;

    XMLElement *pRoot = xmlDoc.NewElement("VTKFile");
    pRoot->SetAttribute("type", "PolyData");
    pRoot->SetAttribute("version", "0.1");
    xmlDoc.InsertFirstChild(pRoot);

    {
        // "Polydata" Sub Tree
        XMLElement *pPolyData = xmlDoc.NewElement("PolyData");

        {
            // "Piece" Sub Tree
            XMLElement *pPiece = xmlDoc.NewElement("Piece");
            // "NumberOfPoints" attribute, indicating the number of points (particles) to draw
            pPiece->SetAttribute("NumberOfPoints", nPart); // Particles in the simulation
                        pPiece->SetAttribute("NumberOfLines", nEdges); // Edges in the simulation

            {
                // Points(Nanoparticles) coordinates
                XMLElement *pPoints = xmlDoc.NewElement("Points");

                {
                    XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
                    pDataArray->SetAttribute("NumberOfComponents", 3); // 3 Dimensions
                    pDataArray->SetAttribute("type", "Float32");
                    pDataArray->SetAttribute("format", "ascii");

                                        std::string s_coords = V_To_S(coords);
                    pDataArray->SetText(s_coords.c_str());


                    pPoints->InsertEndChild(pDataArray);
                }

                pPiece->InsertEndChild(pPoints);
            }

            {
                XMLElement *pPointData = xmlDoc.NewElement("PointData");
                // Mass scalar for each point (nanoparticle)
                pPointData->SetAttribute("Scalars", "particles data");

    {
          // Masses
        XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
        pDataArray->SetAttribute("type", "Float32");
        pDataArray->SetAttribute("Name", "mass");
        pDataArray->SetAttribute("format", "ascii");

                    std::string s_masses = V_To_S(masses);
                    pDataArray->SetText(s_masses.c_str());

                    pPointData->InsertEndChild(pDataArray);
                }

    {
          // particles diameter
        XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
        pDataArray->SetAttribute("type", "Float32");
        pDataArray->SetAttribute("Name", "particles diameter");
        pDataArray->SetAttribute("format", "ascii");

                    std::string s_diameter = V_To_S(diameter);
                    pDataArray->SetText(s_diameter.c_str());

                    pPointData->InsertEndChild(pDataArray);
                }

        {
          // aggregates diameter
          XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
          pDataArray->SetAttribute("type", "Float32");
          pDataArray->SetAttribute("Name", "aggregates diameter");
          pDataArray->SetAttribute("format", "ascii");

          std::string s_agg_diameter = V_To_S(agg_diameter);
          pDataArray->SetText(s_agg_diameter.c_str());

          pPointData->InsertEndChild(pDataArray);
        }

    {
          // aggregate ID in which are included
        XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
        pDataArray->SetAttribute("type", "Int32");
        pDataArray->SetAttribute("Name", "aggregateId");
        pDataArray->SetAttribute("format", "ascii");

                    std::string s_affil = V_To_S(affil);
                    pDataArray->SetText(s_affil.c_str());

                    pPointData->InsertEndChild(pDataArray);
                }

        {
          // particle ID
          XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
          pDataArray->SetAttribute("type", "Int32");
          pDataArray->SetAttribute("Name", "particleId");
          pDataArray->SetAttribute("format", "ascii");

          std::string s_p_IDs = V_To_S(p_IDs);
          pDataArray->SetText(s_p_IDs.c_str());

          pPointData->InsertEndChild(pDataArray);
        }

                pPiece->InsertEndChild(pPointData);
            }

      // Edges connecting the particles
    ////////////////////////////////////////////////////////////////////////
      XMLElement *pCellData = xmlDoc.NewElement("CellData");
      {
        // DataArray "edge_type" indicating the type of the edge, 0 physical, 1 contraints
        {
          XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
          pDataArray->SetAttribute("type", "Int32");
          pDataArray->SetAttribute("Name", "edge_type");
          pDataArray->SetAttribute("format", "ascii");

          std::string s_types = V_To_S(abs_idx[1]);
          pDataArray->SetText(s_types.c_str());

          pCellData->InsertEndChild(pDataArray);
        }

        pPiece->InsertEndChild(pCellData);
      }

      XMLElement *pLines = xmlDoc.NewElement("Lines");
      {
        // DataArray "Connectivity", points indices couples for indicating the edge extremes
        {
          XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
          pDataArray->SetAttribute("type", "Int32");
          pDataArray->SetAttribute("Name", "connectivity");
          pDataArray->SetAttribute("format", "ascii");

          std::string s_edges = V_To_S(abs_idx[0]);
          pDataArray->SetText(s_edges.c_str());

          pLines->InsertEndChild(pDataArray);
        }
        // DataArray "offset"
        {
          XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
          pDataArray->SetAttribute("type", "Float32");
          pDataArray->SetAttribute("Name", "offsets");
          pDataArray->SetAttribute("format", "ascii");

          std::string s_types = V_To_S(abs_idx[2]);
          pDataArray->SetText(s_types.c_str());

          pLines->InsertEndChild(pDataArray);
        }
        pPiece->InsertEndChild(pLines);
      }


                //////////////////////////////////////////////////////////////////////////////////////////
            pPolyData->InsertEndChild(pPiece);
        }

        pRoot->InsertEndChild(pPolyData);
    }


        std::string file_name = root_path + "nanodome" + std::to_string(iter) + ".vtp";
    xmlDoc.SaveFile(file_name.c_str());



}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::motion(double dt, double T) {

  // Initialize the model if not done before
  if (init == false) { initialize(); }

    double cube_side = pow(this->volume, 1. / 3.);

    // For each aggregate in the simulation
    for (auto& agg : this->aggregates) {

        // update the particle positions
        //agg.vsv_update_x_lgv(p.dt, gamma, p.T);
        agg->langevin_x_update(dt,gp,T);

        // update the velocity
        //agg.vsv_update_v_lgv(p.dt);
        agg->langevin_v_update(dt);

        // check reflections on the wall
        agg->wall_reflection(cube_side);
    }

    // apply aggregate constraints
    for (auto& a : this->aggregates)
        a->SHAKE();
}



template<typename A>
double ConstrainedLangevinParticlePhase<A>::nucleation(double dt, double J, double j) {

  // Initialize the model if not done before
  if (init == false) { initialize(); }

    // average number of nucleation events in this timestep
    double n_events = J * this->volume * dt;

    // single nucleation events
    int n = int(n_events);

    // fractional part
    double n_frac = n_events - n;

    // add all the new particles according to J
    for (int i = 0; i<n; ++i)
        this->add_particle(j);

    // the fraction < 1 is treated in a probabilistic way
    if (ndm::uniform_double_distr(ndm::rand_gen) <= n_frac) {
        this->add_particle(j);
        n++;
    }

    // molecules consuption per unit volume and time
    return n*j / (this->volume*dt);
}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::aggregates_number_balance() {

  // Initialize the model if not done before
  if (init == false) { initialize(); }

  int N = this->aggregates.size();

  if (N>this->max_aggregates_number) {

    auto it = this->aggregates.begin();
    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    std::advance(it, int(rho*(N - 1)));

    // ID erased aggregate
    int eraased_ID = (*it)->get_id();
    // center of mass
    std::valarray<double> com = (*it)->get_center_of_mass();

    this->aggregates.erase(it);

    // Erase from the grid
    // Get the cell coordinates
    std::valarray<int> cell_coords = this->grid.get_cell_index(com);
    // erase from the grid
    auto aggregates_cell = this->grid(cell_coords[0], cell_coords[1], cell_coords[2]).get_list();
    int agg_idx = 0;
    for (auto agg = aggregates_cell.begin(); agg != aggregates_cell.end(); agg++) {
      if ((*agg)->get_id() == eraased_ID) {
        this->grid(cell_coords[0], cell_coords[1], cell_coords[2]).erase_n(agg_idx);
      }
      else
        agg_idx++;
    }

    // Update volume
    this->volume *= (N - 1.0) / N;
  }

  if (N < this->min_aggregates_number && N > 0 ) {

    while (N != this->max_aggregates_number) {

      double side = pow(this->volume, 1.0 / 3.0);

      auto it = this->aggregates.begin();
      double rho = ndm::uniform_double_distr(ndm::rand_gen);

      std::advance(it, int(rho*(N - 1)));

      // Using the copy contructor, create a shared pointer for the aggregate copy
      auto agg = std::make_shared<RATTLEAggregate<DynamicParticle>>(*(*it));

      // Debug
      std::cout << "ID Aggregate to copy: "
        << (*it)->get_id()
        << "ID Copied Aggregate: "
        << agg->get_id() << std::endl;

      // get the enclosing sphere radius
      double r0 = agg->get_enclosing_sphere_diameter()*0.5;

      // get volume
      double v0 = agg->get_volume();

      // while guard
      bool placement = false;
      // new center of mass coordinates
      std::valarray<double> n_c_mass;

      while (!placement) {
        // get a random new position to place the aggregate
        n_c_mass = {
          ndm::uniform_double_distr(ndm::rand_gen),
          ndm::uniform_double_distr(ndm::rand_gen),
          ndm::uniform_double_distr(ndm::rand_gen) };

        //center of mass scaled in the domain
        n_c_mass = side*(n_c_mass - 0.5);

        placement = true;

        // Check if in the new position, can collide or is superposed with some aggregate

        // get the cell where the aggregate can be placed
        std::valarray<int> cell_idx = this->grid.get_cell_index(n_c_mass);

        std::vector<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> agg_cell_list;

        agg_cell_list = this->grid(cell_idx[0], cell_idx[1], cell_idx[2]).get_list();

        //try {
        //	// Get the aggregates in the cell
        //	agg_cell_list = this->grid(cell_idx[0], cell_idx[1], cell_idx[2]).get_list();
        //}
        //catch(std::exception& e) {
        //	std::cout << n_c_mass[0] << " " << n_c_mass[1] << " " << n_c_mass[2] << std::endl;
        //	std::cout << cell_idx[0] << " " << cell_idx[1] << " " << cell_idx[2] << std::endl;
        //}

        for (auto a1 = agg_cell_list.begin(); a1 != agg_cell_list.end(); ++a1) {
          double r1 = (*a1)->get_enclosing_sphere_diameter()*0.5;
          std::valarray<double> diff = n_c_mass - (*a1)->get_center_of_mass();
          double dist = sqrt((diff*diff).sum());
          // if the collision is possible, new center
          if (dist - r0 - r1 <= 1.0e-10) {
            std::cout << "POSSIBLE COLLISION IN DUPLICATING AGGREGATE, NEW CENTER OF MASS" << std::endl;
            std::cout << n_c_mass[0] << " " << n_c_mass[1] << " " << n_c_mass[2]<< std::endl;
            std::cout << "Volume: " << v0 << std::endl;
            std::cout << cell_idx[0] << " " << cell_idx[1] << " " << cell_idx[2]
          << " cell side: " << this->grid.get_cell_side() << " volume side: " << side << " grid side: " << this->grid.get_grid_side()
          << std::endl << "Volume_side - grid side = " << abs(side - this->grid.get_grid_side())
          << "Cells: " << this->grid.get_C() << std::endl;

            placement = false;
            break;
          }
        }

      }

      // Change the aggregate's center of mass
      agg->change_center_of_mass(n_c_mass);

      // Add the the aggregates
      this->aggregates.push_back(agg);

      // Insert in the grid
      this->grid.insert_aggregate(agg);

      // Update volume
      this->volume *= (N + 1.0) / N;

      // Update the grid
      int biggest_agg = -1;
      double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
      this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

      N++;
    }

  }

}


template<typename A>
void ConstrainedLangevinParticlePhase<A>::coagulation() {

  // Initialize the model if not done before
  if (this->init == false) { this->initialize(); }

    for (auto a0 = this->aggregates.begin(); a0 != this->aggregates.end(); ++a0) {
        for (auto a1 = std::next(a0); a1 != this->aggregates.end(); ++a1) {
            std::shared_ptr<DynamicParticle> p0;
            std::shared_ptr<DynamicParticle> p1;
            bool collision = (*a0)->check_collision((*a1), 1.0e-10, p0, p1);
            if (collision) {
                // Merge aggregates
                (*a0)->merge((*a1), p0, p1);
                // Update actual particles cardinality inside the agregate (for sintering)
                (*a0)->update_particles_cardinality();
                // Update graph of the new aggregate
                (*a0)->update_graph();
                // Erase graph of the included aggregate
                (*a1)->erase_graph();
#ifdef VERBOSE
                std::cout << "CHECK INDICES:" << std::endl;
                std::cout << "ID 0: " << (*a0)->get_id() << " ID 1: " << (*a1)->get_id() << std::endl;
#endif

                // trick for deleting the current bond without crash
                if (a1 == this->aggregates.begin()) {
                    this->aggregates.remove((*a1));
                    break;
                }
                else {
                    a1 = std::prev(a1);
                    this->aggregates.remove((*(std::next(a1))));
                }

            }


        }
    }

    return;
}


#endif // CONSTRAINEDLANGEVINPARTICLEPHASE_H
