#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <mechanism.H>
#include <PelePhysics.H>
#include <AMReX_GpuDevice.H>

#include <ReactorBase.H>

AMREX_GPU_DEVICE
inline void
initialize_data(
  int i,
  int j,
  int k,
  int testcase,
  amrex::Array4<amrex::Real> const& X,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rhoX
) noexcept
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Xtemp;
  amrex::Real density_cgs;

  amrex::Real pressure = 1013250;
  for (int n = 0; n < NUM_SPECIES; n++){
    if (testcase != 0) X(i, j, k, n) = Xtemp[n] = 0;
    else X(i, j, k, n) = Xtemp[n] = 1/amrex::Real(NUM_SPECIES);
  }

  if (testcase != 0){
     temp(i, j, k) = 900;
     X(i, j, k, CH4_ID) = Xtemp[CH4_ID] = 0.5;
     X(i, j, k,  O2_ID) = Xtemp[ O2_ID] = 0.5;
  }
  else temp(i, j, k) = 300;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PXT2R(pressure, &Xtemp[0], temp(i,j,k), density_cgs);
  rhoX(i,j,k) = density_cgs;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {

    amrex::ParmParse pp;
    int testcase = 0;
    int nc = 1;

    pp.query("testcase", testcase);

    amrex::Box domain(
      amrex::IntVect(D_DECL(0, 0, 0)),
      amrex::IntVect(D_DECL(nc-1, 0, 0)));

    amrex::GpuArray<amrex::Real, 3> plo, phi, dx;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
      phi[i] = domain.length(i);
      dx[i] = (phi[i] - plo[i]) / domain.length(i);
    }

    amrex::RealBox real_box(
      {AMREX_D_DECL(plo[0], plo[1], plo[2])},
      {AMREX_D_DECL(phi[0], phi[1], phi[2])});

    int coord = 0;

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord, is_periodic);

    int max_size = 128;
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    const int num_spec = NUM_SPECIES;

    amrex::DistributionMapping dm{ba};

    int num_grow = 0;
    amrex::MultiFab mole_frac(ba, dm, num_spec, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab molar_density(ba, dm, 1, num_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& gbox = mfi.tilebox();

      amrex::Array4<amrex::Real> const& X_a = mole_frac.array(mfi);
      amrex::Array4<amrex::Real> const& T_a = temperature.array(mfi);
      amrex::Array4<amrex::Real> const& rho_a = molar_density.array(mfi);

      amrex::ParallelFor(
        gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(i, j, k, testcase, X_a, T_a, rho_a);
        });
    }

    amrex::Vector<std::string> species_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(species_names);
    amrex::MultiFab wdots(ba, dm, num_spec, num_grow);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& box = mfi.tilebox();

      const auto X = mole_frac.array(mfi);
      const auto temp = temperature.array(mfi);
      const auto rho = molar_density.array(mfi);
      const auto wdot = wdots.array(mfi);

      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real Xl[NUM_SPECIES];
        amrex::Real Wl[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n){
          Xl[n] = X(i, j, k, n);}
        auto eos = pele::physics::PhysicsType::eos();
        eos.RTX2WDOT( rho(i,j,k), temp(i,j,k), &Xl[0], &Wl[0]);
        for (int n = 0; n < NUM_SPECIES; ++n){
          wdot(i,j,k,n) = Wl[n];
          amrex::Print() << "wdot_" << species_names[n] << " = " << wdot(i,j,k,n) << std::endl;
          }
      }); 
    }
  }
  amrex::Finalize();
  return 0;
}
