#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>

void createGeoemtry(amrex::Box& computationalDomain,
                    amrex::RealBox& physicalDomain,
                    amrex::Geometry& geometry)
{
    // upper and lower bounds of computational domain
    amrex::IntVect lb(AMREX_D_DECL(0,0,0));
    amrex::IntVect ub(AMREX_D_DECL(127,0,0));
    // computational domain
    computationalDomain = amrex::Box(lb, ub);
    // physical domain
    physicalDomain = amrex::RealBox({AMREX_D_DECL(0.0, 0.0, 0.0)},
                                    {AMREX_D_DECL(1.0, 0.0, 0.0)});
    // Boundary conditions in three directions
    amrex::Array<int, AMREX_SPACEDIM> bcs{AMREX_D_DECL(0, 1, 1)};
    geometry.define(computationalDomain, physicalDomain,
                    amrex::CoordSys::cartesian, bcs);
}

int run(){
/****************** Initialize Geometry *******************/
    amrex::Box domain;
    amrex::RealBox realdomain;
    amrex::Geometry geometry;
    createGeoemtry(domain, realdomain, geometry);

    /******************** Initialize blocks *******************/
    int max_grid_size = 128;
    amrex::BoxArray boxarray;
    boxarray.define(domain);
    boxarray.maxSize(max_grid_size);
    /******************** Distribute blocks *******************/
    amrex::DistributionMapping distmap{boxarray};
    /*********** Isolate Spatial Discretization ***************/
    amrex::Real dx = geometry.CellSize(0);
    /***************** Ghost cells around MultiFabs ***********/
    int n_ghost_cells = 0;
    /************* Face Centered MultiFab Arrays **************/
    amrex::MultiFab mole_frac(amrex::convert(
                                    boxarray,
                                    amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, NUM_SPECIES+NUM_SURFACE_SPECIES, n_ghost_cells);
    amrex::MultiFab mass_frac(amrex::convert(
                                    boxarray,
                                    amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, NUM_SPECIES, n_ghost_cells);
    amrex::MultiFab temperature(amrex::convert(
                                    boxarray,
                                    amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    amrex::MultiFab density(amrex::convert(
                                    boxarray,
                                    amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    amrex::MultiFab cp(amrex::convert(
                                boxarray,
                                amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    amrex::MultiFab cv(amrex::convert(
                                boxarray,
                                amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    amrex::MultiFab cp_X(amrex::convert(
                                boxarray,
                                amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    amrex::MultiFab cv_X(amrex::convert(
                                boxarray,
                                amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 1, n_ghost_cells);
    /******************* Set Initial State *******************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& bx = mfi.tilebox();
      auto const& X_a = mole_frac.array(mfi);
      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      amrex::ParallelFor(
        bx, [X_a, Y_a, T_a, dx] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initdata(i, j, k, X_a, Y_a, T_a, dx);
        });
    }
    /*************** Compute gas density given PTY ***************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& rho_a = density.array(mfi);
      amrex::ParallelFor(
        bx, [rho_a, T_a, Y_a] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_density(i, j, k, rho_a, T_a, Y_a);
        });
    }
    /******************* Compute gas Cp given TY ******************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& cp_a = cp.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_cp(i, j, k, Y_a, T_a, cp_a);
        });
    }
    /******************* Compute gas Cv given TY ******************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& cv_a = cv.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_cv(i, j, k, Y_a, T_a, cv_a);
        });
    }
    /*************** Compute surface Cp given TX **************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& X_a = mole_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& cp_a = cp_X.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_surface_cp(i, j, k, X_a, T_a, cp_a);
        });
    }

    /*************** Compute surface Cv given TX **************/
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& X_a = mole_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& cv_a = cv_X.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_surface_cv(i, j, k, X_a, T_a, cv_a);
        });
    }
    /********************** Create PLotfile *********************/

    // Gas phase mass fractions, density, cp, cv
    amrex::MultiFab VarPlt(amrex::convert(
                                boxarray,
                                amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                                distmap, 3, n_ghost_cells);

    amrex::MultiFab::Copy(VarPlt, density, 0, 0, 1, n_ghost_cells);
    amrex::MultiFab::Copy(VarPlt, cv, 0, 1, 1, n_ghost_cells);
    amrex::MultiFab::Copy(VarPlt, cp, 0, 2, 1, n_ghost_cells);

    std::string pltfile("gas");
    std::string outfile = amrex::Concatenate(pltfile, 1);
    amrex::Vector<std::string> plt_VarsName;
    plt_VarsName.push_back("density");
    plt_VarsName.push_back("cv");
    plt_VarsName.push_back("cp");

    amrex::WriteSingleLevelPlotfile(
      outfile, VarPlt, plt_VarsName, geometry, 0.0, 0);

  // Surface cp, cv based on Mole fractions
  amrex::MultiFab surf_VarPlt(amrex::convert(
                              boxarray,
                              amrex::IntVect{AMREX_D_DECL(1,0,0)}),
                              distmap, 2, n_ghost_cells);

  amrex::MultiFab::Copy(surf_VarPlt, cv_X, 0, 0, 1, n_ghost_cells);
  amrex::MultiFab::Copy(surf_VarPlt, cp_X, 0, 1, 1, n_ghost_cells);

  std::string surf_pltfile("surface");
  std::string surf_outfile = amrex::Concatenate(surf_pltfile, 1);
  amrex::Vector<std::string> surf_plt_VarsName;
  surf_plt_VarsName.push_back("cv");
  surf_plt_VarsName.push_back("cp");

  amrex::Vector<std::string> surf_species_names;
  pele::physics::surface::speciesNames<pele::physics::PhysicsType::surface_type>(surf_species_names);

  amrex::WriteSingleLevelPlotfile(
    surf_outfile, surf_VarPlt, surf_plt_VarsName, geometry, 0.0, 0);

  amrex::Print() << "\nSurface Species in Reaction Mechanism:" << std::endl;

  for(int k = 0; k < NUM_SURFACE_SPECIES; ++k)
    amrex::Print() << surf_species_names[k] << "\t";

  amrex::Print() << std::endl;
    return 0;
}


int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  int flag = run();
  amrex::Finalize();
  return 0;
}
