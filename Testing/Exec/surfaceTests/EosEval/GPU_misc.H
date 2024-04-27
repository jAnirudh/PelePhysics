#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

#if defined(NUM_SURFACE_SPECIES)
#define ALL_SPECIES NUM_SPECIES+NUM_SURFACE_SPECIES
#else 
#define ALL_SPECIES NUM_SPECIES
#endif

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(int i, int j, int k,
         amrex::Array4<amrex::Real> const& Y,
         amrex::Array4<amrex::Real> const& T,
         amrex::Real const& dx
        )
{
    amrex::Real pressure = 1013250.0; // dyne/cm^2
    amrex::Real Tinlet = 300.0;       // Kelvin
    amrex::Real Tsurface = 900.0;     // Kelvin
    
    // Initialize linearly varying Temperature profile
    T(i,j,k) = Tinlet + i*dx*(Tsurface-Tinlet);

    // Initialize species
    amrex::GpuArray<amrex::Real, NUM_SPECIES> X;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Ytemp;

    for (int n=0; n < NUM_SPECIES; ++n){
        X[n] = 0;}

    X[H2_ID]  = 0.05;
    X[O2_ID]  = 0.21;
    X[N2_ID]  = 0.78;
    X[AR_ID]  = 0.01;
    X[CH4_ID] = 0.095;

    // normalize gas phase mole-fractions to equal 1
    amrex::Real sum = 0.0;
    for (int n=0; n < NUM_SPECIES; ++n){
        sum += X[n];}
    for (int n=0; n < NUM_SPECIES; ++n){
        X[n] /= sum;}

    auto eos = pele::physics::PhysicsType::eos();
    eos.X2Y(&X[0], &Ytemp[0]);

    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y(i,j,k,n) = Ytemp[n];}
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_density(
    int i,
    int j,
    int k,
    amrex::Array4<amrex::Real> const& rho,
    amrex::Array4<amrex::Real> const& T,
    amrex::Array4<amrex::Real> const& Y) noexcept
{
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Ytemp;
    for (int n = 0; n < NUM_SPECIES; ++n) {
        Ytemp[n] = Y(i, j, k, n);
    }

    amrex::Real pressure = 1013250.0; // dyne/cm^2
    auto eos = pele::physics::PhysicsType::eos();
    // Initialize density
    amrex::Real rho_cgs;
    eos.PYT2R(pressure, &Ytemp[0], T(i,j,k), rho_cgs);
    rho(i,j,k) = rho_cgs * 1.0e3; // CGS -> MKS conversion
}

 
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cp(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& cp) noexcept
{

  amrex::GpuArray<amrex::Real, ALL_SPECIES> mf_pt;
  for (int n = 0; n < ALL_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cp(temp(i, j, k), mf_pt.arr, cp(i, j, k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_cv(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& cv) noexcept
{

  amrex::GpuArray<amrex::Real, ALL_SPECIES> mf_pt;
  for (int n = 0; n < ALL_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cv(temp(i, j, k), mf_pt.arr, cv(i, j, k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_T_from_EY(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<const amrex::Real> const& energy) noexcept
{

  amrex::GpuArray<amrex::Real, ALL_SPECIES> mf_pt;
  for (int n = 0; n < ALL_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.EY2T(energy(i, j, k), mf_pt.arr, temp(i, j, k));
}

#endif
