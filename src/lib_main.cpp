#include "ontodome.h"

//// PyBind11 tools
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
namespace py = pybind11;



/// Internal tools
nanoNetwork create_network(std::vector<std::string> specs, std::vector<double> pp, std::vector<double> TT, std::vector<std::valarray<double>> cs, double mass, double bulk_l) {

  std::vector<nanoCell*> cells;

  for (int k=0; k < int(cs.size()); k++) {
    cells.push_back(new nanoCell(specs, pp[k], TT[k], cs[k], mass, bulk_l));
  }

  nanoNetwork net(cells);

  return net;
}

/// Trampoline methods. Avoids PyBind11 static cast.
template <class R, class O1, class O2>
void create_relation_to(O1* o1, O2* o2) {
  o1->template createRelationTo<R, O2>(o2);
}

template <class O1, class O2>
std::vector<O2*> get_related_objects(O1* o1, O2 o2) {
  return o1->template getRelatedObjects<O2>();
}

template <class O1>
std::vector<std::string> list_relations(O1* o1) {
  std::vector<std::string> res;
  for (auto x : o1->relations) {
      if (x->getRange()->getUuid() != o1->getUuid()) {
          res.push_back(x->getRange()->getClassName());
      } else {
          res.push_back(x->getDomain()->getClassName());
      }
  }
  return res;
}

// Binding module
PYBIND11_MODULE(libontodome, m) {
  m.doc() = "NanoDOME python plugin library for SimPhony";

  py::class_<Real>(m,"Real")
      .def(py::init<double>())
      .def_readwrite("value", &Real::value); // to check if python assignments are correctly updated between object via pointers

  py::class_<Unit>(m,"Unit")
      .def(py::init<std::string>());

  py::class_<Time>(m,"Time")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &Time::get_data, &Time::set_data)
      .def_property("unit", &Time::get_unit, &Time::set_unit);

  py::class_<Mass>(m,"Mass")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &Mass::get_data, &Mass::set_data)
      .def_property("unit", &Mass::get_unit, &Mass::set_unit);

  py::class_<BulkDensityLiquid>(m,"BulkDensityLiquid")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &BulkDensityLiquid::get_data, &BulkDensityLiquid::set_data)
      .def_property("unit", &BulkDensityLiquid::get_unit, &BulkDensityLiquid::set_unit);

  py::class_<BulkDensitySolid>(m,"BulkDensitySolid")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &BulkDensitySolid::get_data, &BulkDensitySolid::set_data)
      .def_property("unit", &BulkDensitySolid::get_unit, &BulkDensitySolid::set_unit);

  py::class_<MeltingPoint>(m,"MeltingPoint")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &MeltingPoint::get_data, &MeltingPoint::set_data)
      .def_property("unit", &MeltingPoint::get_unit, &MeltingPoint::set_unit);

  py::class_<Temperature>(m,"Temperature")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &Temperature::get_data, &Temperature::set_data)
      .def_property("unit", &Temperature::get_unit, &Temperature::set_unit);

  py::class_<TemperatureTimeDerivative>(m,"TemperatureTimeDerivative")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &TemperatureTimeDerivative::get_data, &TemperatureTimeDerivative::set_data)
      .def_property("unit", &TemperatureTimeDerivative::get_unit, &TemperatureTimeDerivative::set_unit);

  py::class_<Pressure>(m,"Pressure")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &Pressure::get_data, &Pressure::set_data)
      .def_property("unit", &Pressure::get_unit, &Pressure::set_unit);

  py::class_<PressureTimeDerivative>(m,"PressureTimeDerivative")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &Pressure::get_data, &Pressure::set_data)
      .def_property("unit", &Pressure::get_unit, &Pressure::set_unit);

  py::class_<SingleComponentComposition>(m,"SingleComponentComposition")
      .def(py::init<MolarFraction*,std::string>())
      .def_readwrite("name", &SingleComponentComposition::name)
      .def("list_relations",list_relations<SingleComponentComposition>)
      .def("create_relation_to",create_relation_to<hasProperty,SingleComponentComposition,SurfaceTension>)
      .def("create_relation_to",create_relation_to<hasProperty,SingleComponentComposition,SaturationPressure>)
      .def("create_relation_to",create_relation_to<hasScalarProperty,SingleComponentComposition,Mass>)
      .def("create_relation_to",create_relation_to<hasScalarProperty,SingleComponentComposition,BulkDensityLiquid>)
      .def("create_relation_to",create_relation_to<hasScalarProperty,SingleComponentComposition,BulkDensitySolid>)
      .def("create_relation_to",create_relation_to<hasScalarProperty,SingleComponentComposition,MeltingPoint>)
      .def("create_relation_to",create_relation_to<hasScalarProperty,SingleComponentComposition,MolarFraction>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,SurfaceTension>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,SaturationPressure>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,Mass>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,BulkDensityLiquid>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,BulkDensitySolid>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,MeltingPoint>)
      .def("get_related_objects",get_related_objects<SingleComponentComposition,MolarFraction>);

  py::class_<GasMixture>(m,"GasMixture")
      .def(py::init<>())
      .def("list_relations",list_relations<GasMixture>)
      .def("create_relation_to",create_relation_to<hasPart,GasMixture,SingleComponentComposition>)
      .def("create_relation_to",create_relation_to<hasProperty,GasMixture,Time>)
      .def("create_relation_to",create_relation_to<hasProperty,GasMixture,Temperature>)
      .def("create_relation_to",create_relation_to<hasProperty,GasMixture,TemperatureTimeDerivative>)
      .def("create_relation_to",create_relation_to<hasProperty,GasMixture,Pressure>)
      .def("create_relation_to",create_relation_to<hasProperty,GasMixture,PressureTimeDerivative>)
      .def("get_related_objects",get_related_objects<GasMixture,SingleComponentComposition>)
      .def("get_related_objects",get_related_objects<GasMixture,Time>)
      .def("get_related_objects",get_related_objects<GasMixture,Temperature>)
      .def("get_related_objects",get_related_objects<GasMixture,TemperatureTimeDerivative>)
      .def("get_related_objects",get_related_objects<GasMixture,Pressure>)
      .def("get_related_objects",get_related_objects<GasMixture,PressureTimeDerivative>);

  py::class_<MolarFraction>(m,"MolarFraction")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &MolarFraction::get_data, &MolarFraction::set_data)
      .def_property("unit", &MolarFraction::get_unit, &MolarFraction::set_unit);

  py::class_<SurfaceTensionPolynomialSoftwareModel>(m,"SurfaceTensionPolynomialSoftwareModel")
      .def(py::init<>());

  py::class_<SurfaceTensionMaterialRelation>(m,"SurfaceTensionMaterialRelation")
      .def(py::init<>())
      .def("run",&SurfaceTensionMaterialRelation::run)
      .def("create_relation_to",create_relation_to<hasSoftwareModel,SurfaceTensionMaterialRelation,SurfaceTensionPolynomialSoftwareModel>);

  py::class_<SurfaceTension>(m,"SurfaceTension")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &SurfaceTension::get_data, &SurfaceTension::set_data)
      .def_property("unit", &SurfaceTension::get_unit, &SurfaceTension::set_unit)
      .def("create_relation_to",create_relation_to<hasMathematicalModel,SurfaceTension,SurfaceTensionMaterialRelation>)
      .def("get_related_objects",get_related_objects<SurfaceTension,SurfaceTensionMaterialRelation>);

  py::class_<SaturationPressurePolynomialSoftwareModel>(m,"SaturationPressurePolynomialSoftwareModel")
      .def(py::init<>());

  py::class_<SaturationPressureMaterialRelation>(m,"SaturationPressureMaterialRelation")
      .def(py::init<>())
      .def("run",&SaturationPressureMaterialRelation::run)
      .def("create_relation_to",create_relation_to<hasSoftwareModel,SaturationPressureMaterialRelation,SaturationPressurePolynomialSoftwareModel>);

  py::class_<SaturationPressure>(m,"SaturationPressure")
      .def(py::init<Real*,Unit*>())
      .def_property("value", &SaturationPressure::get_data, &SaturationPressure::set_data)
      .def_property("unit", &SaturationPressure::get_unit, &SaturationPressure::set_unit)
      .def("create_relation_to",create_relation_to<hasMathematicalModel,SaturationPressure,SaturationPressureMaterialRelation>)
      .def("get_related_objects",get_related_objects<SaturationPressure,SaturationPressureMaterialRelation>);

  py::class_<GasModel>(m,"GasModel")
      .def(py::init<>())
      .def("list_relations",list_relations<GasModel>)
      .def("create_relation_to",create_relation_to<hasModel,GasModel,GasMixture>)
      .def("get_related_objects",get_related_objects<GasModel,GasMixture>)
      .def("get_gas_flux",&GasModel::get_gas_flux)
      .def("get_n",&GasModel::get_n)
      .def("timestep",&GasModel::timestep)
      .def("print",&GasModel::print);

  py::class_<ClassicalNucleationTheory>(m,"ClassicalNucleationTheory")
      .def(py::init<>())
      .def("list_relations",list_relations<ClassicalNucleationTheory>)
      .def("create_relation_to",create_relation_to<hasModel,ClassicalNucleationTheory,SingleComponentComposition>)
      .def("get_related_objects",get_related_objects<ClassicalNucleationTheory,SingleComponentComposition>)
      .def("stable_cluster_diameter",&ClassicalNucleationTheory::stable_cluster_diameter)
      .def("nucleation_rate",&ClassicalNucleationTheory::nucleation_rate);

  py::class_<MomentModelPratsinis>(m,"MomentModelPratsinis")
      .def(py::init<>())
      .def("list_relations",list_relations<MomentModelPratsinis>)
      .def("timestep",&MomentModelPratsinis::timestep)
      .def("get_mean_diameter",&MomentModelPratsinis::get_mean_diameter)
      .def("get_n_density",&MomentModelPratsinis::get_n_density)
      .def("print_lognormal_val",&MomentModelPratsinis::print_lognormal_val)
      .def("create_relation_to",create_relation_to<hasModel,MomentModelPratsinis,SingleComponentComposition>)
      .def("get_related_objects",get_related_objects<MomentModelPratsinis,SingleComponentComposition>);

  py::class_<PBMFractalParticlePhase<PBMAggregate<Particle>>>(m,"PBMFractalParticlePhase")
      .def(py::init<double,double>(), py::arg("D_f")=1.61, py::arg("volume")=5.0e-16)
      .def("timestep",&PBMFractalParticlePhase<PBMAggregate<Particle>>::timestep)
      .def("create_relation_to",create_relation_to<hasModel,PBMFractalParticlePhase<PBMAggregate<Particle>>,SingleComponentComposition>)
      .def("calc_dt",&PBMFractalParticlePhase<PBMAggregate<Particle>>::calc_dt)
      .def("get_volume",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_volume)
      .def("get_aggregates_mean_spherical_diameter",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_aggregates_mean_spherical_diameter)
      .def("get_mean_particles_number",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_mean_particles_number)
      .def("get_particles_mean_diameter",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_particles_mean_diameter)
      .def("get_particles_sizes",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_particles_sizes)
      .def("get_aggregates_sizes",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_aggregates_sizes)
      .def("get_aggregates_number",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_aggregates_number)
      .def("get_aggregates_density",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_aggregates_density)
      .def("get_mean_sintering_level",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_mean_sintering_level)
      .def("get_mean_fractal_dimension",&PBMFractalParticlePhase<PBMAggregate<Particle>>::get_mean_fractal_dimension)
      .def("volume_expansion",&PBMFractalParticlePhase<PBMAggregate<Particle>>::volume_expansion)
      .def("get_related_objects",get_related_objects<PBMFractalParticlePhase<PBMAggregate<Particle>>,SingleComponentComposition>);

  py::class_<ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>>(m,"ConstrainedLangevinParticlePhase")
      .def(py::init<double>(), py::arg("Vol")=5e-19)
      .def("timestep",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::timestep)
      .def("get_volume",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_volume)
      .def("create_relation_to",create_relation_to<hasModel,ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>,SingleComponentComposition>)
      .def("get_aggregates_mean_spherical_diameter",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_aggregates_mean_spherical_diameter)
      .def("get_mean_particles_number",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_mean_particles_number)
      .def("get_particles_sizes",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_particles_sizes)
      .def("get_particles_mean_diameter",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_particles_mean_diameter)
      .def("get_aggregates_number",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_aggregates_number)
      .def("get_aggregates_sizes",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_aggregates_sizes)
      .def("get_aggregates_density",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_aggregates_density)
      .def("get_mean_sintering_level",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_mean_sintering_level)
      .def("get_mean_fractal_dimension",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_mean_fractal_dimension)
      .def("get_related_objects",get_related_objects<ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>,SingleComponentComposition>)
      .def("volume_expansion",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::volume_expansion)
      .def("get_particles_smallest_diameter",&ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>>::get_particles_smallest_diameter);

  py::class_<nanoCell>(m,"nanoCell")
      .def(py::init<std::vector<std::string>, double, double, std::valarray<double>, double, double>())
      .def("print_aggregates_number",&nanoCell::print_aggregates_number);

  py::class_<nanoNetwork>(m,"nanoNetwork")
      .def(py::init<std::vector<nanoCell*>>())
      .def("print_cells",&nanoNetwork::print_cells)
      .def("get_t",&nanoNetwork::get_t)
      .def("get_dt",&nanoNetwork::get_dt)
      .def("get_cell_particles_diameters",&nanoNetwork::get_cell_particles_diameters)
      .def("get_cell_aggregates_diameters",&nanoNetwork::get_cell_aggregates_diameters)
      .def("get_cell_mean_fractal_dimension",&nanoNetwork::get_cell_mean_fractal_dimension)
      .def("get_cell_pbm_volume",&nanoNetwork::get_cell_pbm_volume)
      .def("timestep",&nanoNetwork::timestep);

}
