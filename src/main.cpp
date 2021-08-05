#include "ontodome.h"

int main()
{
    // Setting a clock to keep track of computational time
    WallClock clock;
    clock.start();

    // GasMixture initial state declaration
    GasMixture gas;

    Time t(new Real(0), new Unit("s"));
    gas.createRelationTo<hasProperty,Time>(&t);

    MolarFraction msi(new Real(0.1), new Unit("#"));
    SingleComponentComposition si(&msi,SiliconSymbol::get_symbol());
    MolarFraction mhe(new Real(0.9), new Unit("#"));
    SingleComponentComposition he(&mhe,HeliumSymbol::get_symbol());

    Temperature T(new Real(3000.), new Unit("K"));
    Pressure p(new Real(101325.), new Unit("Pa"));
    PressureTimeDerivative dpdt(new Real(0.), new Unit("Pa/s"));
    TemperatureTimeDerivative dTdt(new Real(-1e+7), new Unit("K/s"));

    gas.createRelationsTo<hasPart,SingleComponentComposition>({&si,&he});
    gas.createRelationTo<hasProperty,Temperature>(&T);
    gas.createRelationTo<hasProperty,Pressure>(&p);
    gas.createRelationTo<hasProperty,PressureTimeDerivative>(&dpdt);
    gas.createRelationTo<hasProperty,TemperatureTimeDerivative>(&dTdt);

    // Surface Tension and Saturation Pressure settings
    SurfaceTensionPolynomialSoftwareModel stpm;
    SurfaceTensionMaterialRelation stmr;
    SurfaceTension st(new Real(0.), new Unit("N/m"));

    stmr.createRelationTo<hasSoftwareModel,SurfaceTensionPolynomialSoftwareModel>(&stpm);
    st.createRelationTo<hasMathematicalModel,SurfaceTensionMaterialRelation>(&stmr);
    si.createRelationTo<hasProperty,SurfaceTension>(&st);
    si.getRelatedObjects<SurfaceTension>()[0]->getRelatedObjects<SurfaceTensionMaterialRelation>()[0]->run();

    SaturationPressurePolynomialSoftwareModel sapm;
    SaturationPressureMaterialRelation samr;
    SaturationPressure sa(new Real(0.), new Unit("Pa"));

    samr.createRelationTo<hasSoftwareModel,SaturationPressurePolynomialSoftwareModel>(&sapm);
    sa.createRelationTo<hasMathematicalModel,SaturationPressureMaterialRelation>(&samr);
    si.createRelationTo<hasProperty,SaturationPressure>(&sa);

    // Si properties
    Mass sim(new Real(28.0855*AMU), new Unit("kg"));
    sim.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    BulkDensityLiquid sibl(new Real(2570.), new Unit("kg/m3"));
    sibl.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    BulkDensitySolid sibs(new Real(2329.), new Unit("kg/m3"));
    sibs.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    MeltingPoint melp(new Real(1687.), new Unit("K"));
    melp.createRelationTo<hasScalarProperty,SingleComponentComposition>(&si);

    // Helium properties
    Mass hem(new Real(4.005*AMU), new Unit("kg"));
    hem.createRelationTo<hasScalarProperty,SingleComponentComposition>(&he);

    // Models settings
    Time dt(new Real(1e-7), new Unit("s")); // simulation timestep. This definition is not ontologically correct

    GasModel gm;
    gm.createRelationTo<hasModel,GasMixture>(&gas);

    ClassicalNucleationTheory cnt;
    cnt.createRelationTo<hasModel,SingleComponentComposition>(&si);

    MomentModelPratsinis mom;
    mom.createRelationTo<hasModel,SingleComponentComposition>(&si);

//    PBMFractalParticlePhase<PBMAggregate<Particle>> pp(1.61, 5.0e-16);
//    pp.createRelationTo<hasModel,SingleComponentComposition>(&si);

    // Common settings
    int PRINT_EVERY = 1;
    int iter = 0;

    ConstrainedLangevinParticlePhase<RATTLEAggregate<DynamicParticle>> cgmd(5e-19);
    cgmd.createRelationTo<hasModel,SingleComponentComposition>(&si);

    // CGMD loop
    double SAVE_SNAPSHOT = 1.0e-8;
    std::string vtk_path = "E/vtk/";
    double snap_count = 0.0; // Counter for saving VTK file

    // Temporary assignments. To be removed during second CGMD release
    double* T_melt = si.getRelatedObjects<MeltingPoint>()[0]->onData();
    double* s_bulk_density_sol = si.getRelatedObjects<BulkDensitySolid>()[0]->onData();
    double* s_bulk_density_liq = si.getRelatedObjects<BulkDensityLiquid>()[0]->onData();


    // Simulation Main Cycle
    while (*t.onData() <= 0.005) {

        if (gm.get_T() < 300.)
        { *dTdt.onData() = 0.; }

        //Bulk density calculation
        // Initialize the model if not done before
        double bdens = (gm.get_T() < *T_melt) ? *s_bulk_density_sol : *s_bulk_density_liq;

        // check the smallest particle in the system
        double d_min = cgmd.get_particles_smallest_diameter();

        // calculate dt max for langevin equation stability dt (dt < 2*m/alpha)
        double dt_max_lang = d_min*bdens/gm.get_gas_flux();

        // calculate dt max to have v*dt < d/2 for the smallest particle
        double dt_max_coll = sqrt(M_PI*bdens*pow(d_min,5)/(24*3*K_BOL*gm.get_T()));

        if(cgmd.get_aggregates_number()>=1) {
            *dt.onData() = std::min(dt_max_coll,dt_max_lang);
            if(*dt.onData()>1e-10) *dt.onData() = 1.0e-10;
          }

        // particle phase timestep
        double g_cons = cgmd.timestep(*dt.onData(), &gm, &cnt, *T.onData(), &si);

        // gas phase time step
        gm.timestep(*dt.onData(), { -g_cons, 0.0 });

        // Update elapsed time and iterations
        *t.onData() += *dt.onData();
        snap_count += *dt.onData();
        iter++;

        // Print snapshot
        if (counter_trigger(iter,PRINT_EVERY)) {
            std::cout << "t[s]: " << *t.onData() << '\t'
                      << dt_max_coll << '\t'
                      << dt_max_lang << '\t'
                      << "T[K]: " << gm.get_T() << '\t'
                      << "V[m^3]"<< cgmd.get_volume() << '\t'
                      << "|N|: " << cgmd.get_aggregates_number() << '\t'
                      << "max|N|: " << cgmd.get_max_agg_number() << '\t'
                      << "min|N|: " << cgmd.get_min_agg_number() << '\t'
//                      << "|C|: " << cgmd.get_aggregates_cardinality()
                      << std::endl;
          }

//        // Save VTK
//        if (snap_count >= SAVE_SNAPSHOT && pp.get_aggregates_number() > 0) {
//            cgmd.save_vtk(iterations, vtk_path);
//            snap_count = 0.0;
//          }

      }

/*
    // PBM loop
    // loop over timesteps
    while(*t.onData() < 0.005) {

      if (gm.get_T() < 300.)
      { *dTdt.onData() = 0.; }

      // species source term for the gas phase
      double g_si = 0.0;

      // calculate the timestep using an exponential waiting time
      double R_tot = pp.get_total_processes_rate(&gm,&cnt);
      double rho = ndm::uniform_double_distr(ndm::rand_gen);

      // exponential waiting time
      double dt_ = -log(rho)/R_tot;

      // Strang first step
      gm.timestep(dt_/2.0,{0.0,0.0});
      pp.volume_expansion(dt_/2.0,&gm);

      // Strang second step
      g_si += pp.timestep(dt_,&gm,&cnt,&si);
      gm.timestep(dt_,{-g_si,0.0});

      // Strang third step
      gm.timestep(dt_/2.0,{0.0,0.0});
      pp.volume_expansion(dt_/2.0,&gm);

      *t.onData() += dt_;

      iter++;
      if(counter_trigger(iter,PRINT_EVERY)) {

          clock.stop();

          std::cout << "Time: " << *t.onData() << '\t'				// time
                    << "T: " << gm.get_T() << '\t'					// temperature
                    << "S: " << gm.get_S(&si) << '\t'				// supersaturation (S)
                    << "NR: " << cnt.nucleation_rate() << '\t'			// J
                    << "n: " << gm.get_n() << '\t'					// ns
                    << "SCD: " << cnt.stable_cluster_diameter() << '\t'		// j
                    << "MPN: " << pp.get_mean_particles_number() << '\t'		// N_m
                    << "MSL: " << pp.get_mean_sintering_level() << '\t'		//
                    << "AMSD: " << pp.get_aggregates_mean_spherical_diameter() << '\t'
                    << "AN: " << pp.get_aggregates_number() << '\t'
                    << "AD: " << pp.get_aggregates_density() << '\t'
                    << "V: " << pp.get_volume() << '\t'
                    << "MFD: " << pp.get_mean_fractal_dimension() << '\t'
                    << "CI: " << clock.interval()/PRINT_EVERY << std::endl;

          clock.start();
      }
      }

*/
    // Moments loop
/*    while ( *t.onData() <= 0.02) {
      if ( gm.get_T() < 600.) {
        *dTdt.onData() = 0.;
      }

      double g_cons = mom.timestep(*dt.onData());

      gm.timestep(*dt.onData(), {g_cons, 0.0});

      *t.onData() += *dt.onData();
      iter += 1;

      if (counter_trigger(iter, PRINT_EVERY)) {

          auto spec_name = si.name;
          std::cout
              << "Time= "<< *t.onData() << '\t'
              << "Temp= " << gm.get_T() << '\t'
              << "Sat_" << spec_name << "= " << gm.get_S(&si) << '\t'
              << spec_name << "_cons= " << g_cons << '\t'
              << "Mean Diam_" << spec_name << "= " << mom.get_mean_diameter() << '\t'
              << "M0_" << spec_name << "= " << mom.get_n_density() << '\t'
              << "M1_" << spec_name << "= " << mom.get_M1() << '\t'
              << "M2_" << spec_name << "= " << mom.get_M2() << '\t'
              << std::endl << std::endl;
      }
    }
*/
    gm.print();

    clock.stop();
    std::cout << "Execution time: " << clock.interval() << " s" << std::endl;

    return 0;
}
