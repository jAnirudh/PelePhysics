#ifndef _SURFACE_H_
#define _SURFACE_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "PhysicsConstants.H"
#include "IdealSurface.H"

namespace pele {
namespace physics {

#ifdef USE_IDEAL_SURFACE
using SurfaceType = surface::Ideal;
#else
static_assert(false, "Invalid SurfaceType specified");
#endif

namespace surface {

template <typename SurfaceType>
void atomic_weightsCHON(amrex::Real atwCHON[]);

template <typename SurfaceType>
void element_compositionCHON(int ecompCHON[]);

template <typename SurfaceType>
void speciesNames(amrex::Vector<std::string>& spn);

} // namespace SURFACE
} // namespace physics
} // namespace pele
#endif
