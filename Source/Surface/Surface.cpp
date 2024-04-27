#include "Surface.H"

namespace pele {
namespace physics {
namespace surface{

template <typename SurfaceType>
void
speciesNames(amrex::Vector<std::string>& /*spn*/)
{
}

template <>
void
speciesNames<Ideal>(amrex::Vector<std::string>& spn)
{
  SKSYMS_STR(spn);
}

} // namespace surface
} // namespace physics
} // namespace pele
