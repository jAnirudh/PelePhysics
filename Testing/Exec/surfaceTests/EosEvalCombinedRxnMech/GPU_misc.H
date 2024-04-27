#ifndef GPU_MISC_H_
#define GPU_MISC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>
#include <cmath>

#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(int i, int j, int k,
        amrex::Array4<amrex::Real> const& X,
        amrex::Array4<amrex::Real> const& Y,
        amrex::Array4<amrex::Real> const& T,
        amrex::Real const& dx)
{
    amrex::Real pressure = 1013250.0; // dyne/cm^2
    amrex::Real Tinlet = 300.0;       // Kelvin
    amrex::Real Tsurface = 900.0;     // Kelvin
    
    // Initialize linearly varying Temperature profile
    T(i,j,k) = Tinlet + i*dx*(Tsurface-Tinlet);

    // Initialize species

    amrex::GpuArray<amrex::Real, NUM_SPECIES> Xtemp;
    amrex::GpuArray<amrex::Real, NUM_SPECIES> Ytemp;

    for (int n=0; n < NUM_SPECIES; ++n){
        Xtemp[n] = 0; Ytemp[n] = 0;}

    Xtemp[H2_ID]  = 0.05;
    Xtemp[O2_ID]  = 0.21;
    Xtemp[N2_ID]  = 0.78;
    Xtemp[AR_ID]  = 0.01;
    Xtemp[CH4_ID] = 0.095;

    // normalize gas phase mole-fractions to equal 1
    amrex::Real sum = 0.0;
    for (int n=0; n < NUM_SPECIES; ++n){
        sum += Xtemp[n];}
    for (int n=0; n < NUM_SPECIES; ++n){
        Xtemp[n] /= sum;}

    for (int n = 0; n < NUM_SPECIES; ++n) {
        X(i,j,k,n) = Xtemp[n];}

    auto eos = pele::physics::PhysicsType::eos();
    eos.X2Y(&Xtemp[0], &Ytemp[0]);

    for (int n = 0; n < NUM_SPECIES; ++n) {
        Y(i,j,k,n) = Ytemp[n];}

  if (i == 128){
      amrex::GpuArray<amrex::Real, NUM_SURFACE_SPECIES> coverages;

      for (int n=0; n < NUM_SURFACE_SPECIES; ++n)
          coverages[n] = 0;

      coverages[PT_S_ID - NUM_SPECIES] = 6.37367505e-01;
      coverages[H_S_ID - NUM_SPECIES] = 1.07387019e-05;
      coverages[H2O_S_ID - NUM_SPECIES] = 3.82791279e-05;
      coverages[OH_S_ID - NUM_SPECIES] = 1.13297955e-02;
      coverages[CO_S_ID - NUM_SPECIES] = 3.99433519e-02;
      coverages[CO2_S_ID - NUM_SPECIES] = 1.55173855e-07;
      coverages[CH3_S_ID - NUM_SPECIES] = 2.74577790e-07;
      coverages[CH2_Ss_ID - NUM_SPECIES] = 2.74577790e-07;
      coverages[CH_S_ID - NUM_SPECIES] = 2.74577790e-07;
      coverages[C_S_ID - NUM_SPECIES] = 1.71435680e-04;
      coverages[O_S_ID - NUM_SPECIES] = 3.11137915e-01;

    // normalize coverages to 1
    sum = 0.0;
    for (int n=0; n < NUM_SURFACE_SPECIES; ++n){
        sum += coverages[n];}
    for (int n=0; n < NUM_SURFACE_SPECIES; ++n){
        coverages[n] /= sum;}

    for (int n = 0; n < NUM_SURFACE_SPECIES; ++n)
           X(i, j, k, n+NUM_SPECIES) = coverages[n];
   }
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

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
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

  amrex::GpuArray<amrex::Real, NUM_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto eos = pele::physics::PhysicsType::eos();
  eos.TY2Cv(temp(i, j, k), mf_pt.arr, cv(i, j, k));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_surface_cp(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& cp) noexcept
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES+NUM_SURFACE_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES+NUM_SURFACE_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto surface = pele::physics::PhysicsType::surface();
  surface.TX2Cp(temp(i, j, k), mf_pt.arr, cp(i, j, k));

  if (i == 128)
    amrex::Print() << "\nSurface Cp = " << cp(i,j,k) << std::endl;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
get_surface_cv(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& mf,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& cv) noexcept
{

  amrex::GpuArray<amrex::Real, NUM_SPECIES+NUM_SURFACE_SPECIES> mf_pt;
  for (int n = 0; n < NUM_SPECIES+NUM_SURFACE_SPECIES; n++) {
    mf_pt[n] = mf(i, j, k, n);
  }

  auto surface = pele::physics::PhysicsType::surface();
  surface.TX2Cv(temp(i, j, k), mf_pt.arr, cv(i, j, k));

  if (i == 128)
    amrex::Print() << "\nSurface Cv = " << cv(i,j,k) << std::endl;
}

#endif
