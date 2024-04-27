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
  amrex::Array4<amrex::Real> const& X,
  amrex::Array4<amrex::Real> const& cov,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rhoX
) noexcept
{
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Xtemp;
  amrex::Real density_cgs;

  temp(i, j, k) = 900;
  amrex::Real pressure = 1013250;
  for (int n = 0; n < NUM_SPECIES; n++)
    X(i, j, k, n) = Xtemp[n] = 1/amrex::Real(NUM_SPECIES);

  cov(i,j,k,PT_S_ID      - NUM_SPECIES) = 2.53046371e-01;
  cov(i,j,k,H_S_ID       - NUM_SPECIES) = 2.32949269e-02;
  cov(i,j,k,H2O_S_ID     - NUM_SPECIES) = 6.33989436e-05;
  cov(i,j,k,OH_S_ID      - NUM_SPECIES) = 9.42970530e-05;
  cov(i,j,k,CO_S_ID      - NUM_SPECIES) = 6.87250039e-01;
  cov(i,j,k,CO2_S_ID     - NUM_SPECIES) = 1.89501177e-10;
  cov(i,j,k,CH3_S_ID     - NUM_SPECIES) = 1.03735955e-08;
  cov(i,j,k,CH2_Ss_ID    - NUM_SPECIES) = 1.03735955e-08;
  cov(i,j,k,CH_S_ID      - NUM_SPECIES) = 1.03735955e-08;
  cov(i,j,k,C_S_ID       - NUM_SPECIES) = 3.62288519e-02;
  cov(i,j,k,O_S_ID       - NUM_SPECIES) = 2.20839132e-05;

  auto surface = pele::physics::PhysicsType::surface();
  surface.PT2MolarDensity(pressure, temp(i,j,k), density_cgs);
  rhoX(i,j,k) = density_cgs;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {

    amrex::ParmParse pp;
    int nc = 1;
    pp.query("nc", nc);

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
    pp.query("max_size", max_size);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    const int num_spec = NUM_SPECIES;
    const int num_surf_spec = NUM_SURFACE_SPECIES;

    amrex::DistributionMapping dm{ba};

    int num_grow = 0;
    amrex::MultiFab mole_frac(ba, dm, num_spec, num_grow);
    amrex::MultiFab coverages(ba, dm, num_surf_spec, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab molar_density(ba, dm, 1, num_grow);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& gbox = mfi.tilebox();

      amrex::Array4<amrex::Real> const& X_a = mole_frac.array(mfi);
      amrex::Array4<amrex::Real> const& cov_a = coverages.array(mfi);
      amrex::Array4<amrex::Real> const& T_a = temperature.array(mfi);
      amrex::Array4<amrex::Real> const& rho_a = molar_density.array(mfi);

      amrex::ParallelFor(
        gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(i, j, k, X_a, cov_a, T_a, rho_a);
        });
    }

    amrex::Vector<std::string> species_names, surf_species_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(species_names);
    pele::physics::surface::speciesNames<pele::physics::PhysicsType::surface_type>(surf_species_names);

    amrex::MultiFab wdots_gas(ba, dm, num_spec, num_grow);
    amrex::MultiFab wdots_surface(ba, dm, num_spec+num_surf_spec, num_grow);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mole_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& box = mfi.tilebox();

      const auto X = mole_frac.array(mfi);
      const auto cov = coverages.array(mfi);
      const auto temp = temperature.array(mfi);
      const auto rho = molar_density.array(mfi);
      const auto w_g = wdots_gas.array(mfi);
      const auto w_s = wdots_surface.array(mfi);

      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real Xl[num_spec];
        amrex::Real covl[num_surf_spec];
        amrex::Real Wsurface[num_spec + num_surf_spec];
        amrex::Real Wgas[num_spec];
        for (int n = 0; n < num_spec; ++n){
          Xl[n] = X(i, j, k, n);}
        for (int n = 0; n < num_surf_spec; ++n){
          covl[n] = cov(i, j, k, n);}
        auto surface = pele::physics::PhysicsType::surface();
        surface.RTX2WDOTX( rho(i,j,k),
                           temp(i,j,k),
                           &Xl[0],
                           &covl[0],
                           &Wgas[0],
                           &Wsurface[0]);
        for (int n = 0; n < NUM_SPECIES; ++n){
          w_g(i,j,k,n) = Wgas[n];
          w_s(i,j,k,n) = Wsurface[n];
          amrex::Print() << "Gas_wdot_" << species_names[n] << " = "
	                 << w_g(i,j,k,n) << " mol/cm^3/s"
			 << std::endl;
        }
        amrex::Print() << std::endl;
        for (int n = NUM_SPECIES; n < NUM_SPECIES+NUM_SURFACE_SPECIES; ++n)
          w_s(i,j,k,n) = Wsurface[n];
        for (int n = 0; n < NUM_SPECIES+NUM_SURFACE_SPECIES; ++n){
          if (n >= NUM_SPECIES) 
              amrex::Print() << "Surf_wdot_" << surf_species_names[n-NUM_SPECIES]
		             << " = "
			     << w_s(i,j,k,n) << " mol/cm^2/s"
			     << std::endl;
          else
              amrex::Print() << "Surf_wdot_" << species_names[n] << " = "
		             << w_s(i,j,k,n) << " mol/cm^3/s"
			     << std::endl;
        }
      });
    }
  }
  amrex::Finalize();
  return 0;
}
