#include "unitconv.h"
#include <cmath>
#include <boost/format.hpp>

using namespace Pizza;
using namespace std;

const double units::c_SI      = 299792458.0;
const double units::G_SI      = 6.673e-11;
const double units::M_sun_SI  = 1.98892e30;

units units::si()
{
  return units(1.0,1.0,1.0);
}

units units::geom_ulength(double ulength)
{
  return units(ulength, ulength/c_SI, ulength * c_SI * c_SI / G_SI);
}

units units::geom_udensity(double udensity)
{
  return geom_ulength( c_SI / sqrt(G_SI*udensity));
}

units units::geom_umass(double umass)
{
  return geom_ulength(umass*G_SI/(c_SI*c_SI));
}

units units::operator/(const units &base)
{
  return units(ulength/base.ulength, utime/base.utime, umass/base.umass);
}

string units::to_str() const
{
  boost::format fmt("ulength=%.15e m, utime=%.15e s, umass=%.15e kg");
  fmt % ulength % utime % umass;
  return fmt.str();
}

ostream& Pizza::operator<<(ostream& o, const units& u)
{
  o << u.to_str();
  return o;
}


