#ifndef _IDEAL_SURFACE_H_
#define _IDEAL_SURFACE_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include "mechanism.H"
#include "PhysicsConstants.H"

namespace pele {
namespace physics {
namespace surface{

struct Ideal
{
  using surface_type = Ideal;

  static std::string identifier() { return "Ideal"; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void molecular_weight(amrex::Real mw[]) { get_mw(mw); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void inv_molecular_weight(amrex::Real imw[]) { get_imw(imw); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TX2Cp(amrex::Real T, amrex::Real X[], amrex::Real& Cp)
  {
    SKCPBL(&T, X, &Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TX2Cv(amrex::Real T, amrex::Real X[], amrex::Real& Cp)
  {
    SKCVBL(&T, X, &Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
  {
    CKCPBS(&T, Y, &Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2Cp(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
  {
    TY2Cp(T, Y, Cp);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2Cv(amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
  {
    CKCVBS(&T, Y, &Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2Cv(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
  {
    TY2Cv(T, Y, Cv);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Cpi(amrex::Real T, amrex::Real Cpi[]) { CKCPMS(&T, Cpi); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real Y[], amrex::Real& Cs)
  {
    amrex::Real tmp[NUM_SURFACE_SPECIES];
    amrex::Real wbar;
    CKMMWY(&Y[0], &wbar);
    amrex::Real T = P * wbar / (R * Constants::RU);
    CKCVMS(&T, tmp);
    amrex::Real Cv = 0.0;
    for (int i = 0; i < NUM_SURFACE_SPECIES; i++)
      Cv += Y[i] * tmp[i];
    amrex::Real G = (wbar * Cv + Constants::RU) / (wbar * Cv);
    Cs = std::sqrt(G * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
  {
    amrex::Real tmp[NUM_SURFACE_SPECIES];
    amrex::Real P;
    CKPY(&R, &T, Y, &P);
    CKCVMS(&T, tmp);
    amrex::Real Cv = 0.0;
    for (int i = 0; i < NUM_SURFACE_SPECIES; ++i)
      Cv += Y[i] * tmp[i];
    amrex::Real wbar;
    CKMMWY(&Y[0], &wbar);
    amrex::Real G = (wbar * Cv + Constants::RU) / (wbar * Cv);
    Cs = std::sqrt(G * P / R);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
  {
    // For Fuego this function is really just a wrapper for GET_T_GIVEN_EY
    // In SRK this will be different probably
    int lierr = 0;
    GET_T_GIVEN_EY(&E, Y, &T, &lierr);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void REY2T(amrex::Real /*R*/, amrex::Real E, amrex::Real Y[], amrex::Real& T)
  {
    EY2T(E, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void HY2T(amrex::Real H, amrex::Real Y[], amrex::Real& T)
  {
    // For Fuego this function is really just a wrapper for GET_T_GIVEN_HY
    // In SRK this will be different probably
    int lierr = 0;
    GET_T_GIVEN_HY(&H, Y, &T, &lierr);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RHY2T(amrex::Real /*R*/, amrex::Real H, amrex::Real Y[], amrex::Real& T)
  {
    HY2T(H, Y, T);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYET2P(
    amrex::Real R,
    amrex::Real Y[],
    amrex::Real& E,
    amrex::Real& T,
    amrex::Real& P)
  {
    int lierr = 0;
    GET_T_GIVEN_EY(&E, Y, &T, &lierr);
    CKPY(&R, &T, Y, &P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P)
  {
    CKPY(&R, &T, Y, &P);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2T(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& T)
  {
    amrex::Real wbar;
    CKMMWY(&Y[0], &wbar);
    T = P * wbar / (R * Constants::RU);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  RTY2WDOT(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real WDOT[])
  {
    // CKWYR(&R, &T, Y, WDOT);
    amrex::Real C[NUM_SURFACE_SPECIES];
    CKYTCR(&R, &T, Y, C);
    CKWC(&T, C, WDOT);

    amrex::Real mw[NUM_SURFACE_SPECIES];
    get_mw(mw);
    for (int n = 0; n < NUM_SURFACE_SPECIES; n++)
      WDOT[n] *= mw[n];
  }

/** @brief Get molar concentrations given molar gas density,
 *         mole fractions and surface coverages
 *
 *  @param molar_density  molar gas density in mol/cm^3
 *  @param X              array of mole fractions
 *  @param coverages      array of surface coverages
 *  @param concentrations array of molar concentrations
                          (mol/cm^3 for the bulk phase &
                           mol/cm^2 for the interface)
 */
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  getMolarConcentrations(amrex::Real molar_density,
                         amrex::Real X[],
                         amrex::Real coverages[],
                         amrex::Real concentrations[])
  {
   for (int k = 0; k < NUM_SPECIES; ++k)
       concentrations[k] = molar_density * X[k];
   for (int k = 0; k < NUM_SURFACE_SPECIES; ++k)
       concentrations[NUM_SPECIES+k] = SITE_DENSITY * coverages[k];
  }

/** @brief Compute molar production rates given density,
 *         temperature, mole fractions and surface coverages.
 *
 *  @param R      molar gas density in mol/cm^3
 *  @param T      temperature in K
 *  @param X      array of mole fractions
 *  @param COV    array of surface coverages
 *  @param WDOT   array of gas phase molar production rates (mol/cm^3/s)
 *  @param WDOT_S array of molar production rates due to surface reactions
                  (mol/cm^3/s for first NUM_SPECIES and mol/cm^2/s for
                   remaining NUM_SURFACE_SPECIES)
 *  @return void
 */
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  RTX2WDOTX(amrex::Real R,
            amrex::Real T,
            amrex::Real X[],
            amrex::Real COV[],
            amrex::Real WDOT[],
            amrex::Real WDOT_S[])
  {
    // molar concentration arrays
    amrex::Real C[NUM_SPECIES+NUM_SURFACE_SPECIES];
    // compute molar concentrations
    getMolarConcentrations(R, X, COV, C); // from this file
    CKWC(&T, C, WDOT_S, WDOT);
  }

/** @brief Convert molar production rates (mol/cm^3/s) to
 *         production rates (gm/cm^3/s).
 *         WARNING: this method overwrites the original data
 *
 *  @param WDOT   array of gas phase molar production rates (mol/cm^3/s)
 *  @return void
 */
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  MolarWdot2Wdot(amrex::Real WDOT[])
  {
    amrex::Real mw[NUM_SPECIES];
    get_mw(mw);
    for (int n = 0; n < NUM_SPECIES; n++)
      WDOT[n] *= mw[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2JAC(
    amrex::Real R,
    amrex::Real T,
    amrex::Real Y[],
    amrex::Real JAC[],
    const int HP)
  {
    amrex::Real C[NUM_SURFACE_SPECIES];
    CKYTCR(&R, &T, Y, C);
    DWDOT(JAC, C, &T, &HP);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2C(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real acti[])
  {
    CKYTCR(&R, &T, Y, acti);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Ei(amrex::Real T, amrex::Real Ei[]) { CKUMS(&T, Ei); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Ei(
    amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y[]*/, amrex::Real Ei[])
  {
    T2Ei(T, Ei);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2X(amrex::Real Y[], amrex::Real X[]) { CKYTX(Y, X); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void X2Y(amrex::Real X[], amrex::Real Y[]) { CKXTY(X, Y); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2RE(
    amrex::Real P,
    amrex::Real Y[],
    amrex::Real T,
    amrex::Real& R,
    amrex::Real& E)
  {
    CKRHOY(&P, &T, Y, &R);
    amrex::Real ei[NUM_SURFACE_SPECIES];
    CKUMS(&T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SURFACE_SPECIES; n++)
      E += Y[n] * ei[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void PYT2R(amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real& R)
  {
    CKRHOY(&P, &T, Y, &R);
  }

/** @brief Compute molar density given pressure and temperature.
 *
 *  @param P the pressure in dyne/cm^2
 *  @param T the temperature in K
 *  @param R variable address to store density in mol/cm^3
 *  @return void
 */
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void
  PT2MolarDensity(const amrex::Real P,
                  const amrex::Real T,
                  amrex::Real& R)
  {
    R = P / (8.31446261815324e+07 * (T)); // rho = P/(R*T)
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RYP2E(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E)
  {
    amrex::Real wbar;
    CKMMWY(&Y[0], &wbar);
    amrex::Real T = P * wbar / (R * Constants::RU);
    amrex::Real ei[NUM_SURFACE_SPECIES];
    CKUMS(&T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SURFACE_SPECIES; n++)
      E += Y[n] * ei[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2E(amrex::Real T, amrex::Real Y[], amrex::Real& E)
  {
    amrex::Real ei[NUM_SURFACE_SPECIES];
    T2Ei(T, ei);
    E = 0.0;
    for (int n = 0; n < NUM_SURFACE_SPECIES; n++)
      E += Y[n] * ei[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2E(amrex::Real /*R*/, amrex::Real T, amrex::Real Y[], amrex::Real& E)
  {
    TY2E(T, Y, E);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void T2Hi(amrex::Real T, amrex::Real Hi[]) { CKHMS(&T, Hi); }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2Hi(
    amrex::Real /*R*/, amrex::Real T, amrex::Real* /*Y[]*/, amrex::Real Hi[])
  {
    T2Hi(T, Hi);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void S(amrex::Real& s) { s = 1.0; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2G(amrex::Real T, amrex::Real Y[NUM_SURFACE_SPECIES], amrex::Real& G)
  {
    amrex::Real wbar, Cv, Cvx;
    TY2Cv(T, Y, Cv);
    CKMMWY(&Y[0], &wbar);
    Cvx = wbar * Cv;
    G = (Cvx + Constants::RU) / Cvx;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2G(
    amrex::Real /*R*/,
    amrex::Real T,
    amrex::Real Y[NUM_SURFACE_SPECIES],
    amrex::Real& G)
  {
    TY2G(T, Y, G);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void TY2H(amrex::Real T, amrex::Real Y[NUM_SURFACE_SPECIES], amrex::Real& H)
  {
    amrex::Real hi[NUM_SURFACE_SPECIES];
    CKHMS(&T, hi);
    H = 0.0;
    for (int n = 0; n < NUM_SURFACE_SPECIES; n++)
      H += Y[n] * hi[n];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void Y2WBAR(amrex::Real Y[], amrex::Real& WBAR)
  {
    amrex::Real tmp[NUM_SURFACE_SPECIES];
    get_imw(tmp);
    amrex::Real summ = 0.0;
    for (int i = 0; i < NUM_SURFACE_SPECIES; i++)
      summ += Y[i] * tmp[i];
    WBAR = 1.0 / summ;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RPE2dpdr_e(
    amrex::Real R, amrex::Real P, amrex::Real /* E */, amrex::Real& dpdr_e)
  {
    dpdr_e = P / R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RG2dpde(amrex::Real R, amrex::Real G, amrex::Real& dpde)
  {
    dpde = (G - 1.0) * R;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void RTY2dpde_dpdre(
    amrex::Real R,
    amrex::Real T,
    amrex::Real Y[],
    amrex::Real& dpde,
    amrex::Real& dpdr_e)
  {
    amrex::Real G, P;
    TY2G(T, Y, G);
    RTY2P(R, T, Y, P);
    dpde = (G - 1.0) * R;
    dpdr_e = P / R;
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE Ideal(Args...)
  {
  }
};
} // namespace surface
} // namespace physics
} // namespace pele
#endif
