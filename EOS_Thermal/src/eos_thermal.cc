#include "eos_thermal.h"
#include <limits>
#include <stdexcept>
#include <sstream>
#include <iomanip>

using namespace whizza;
using namespace std;

std::ostream& operator<<(std::ostream &s, const eos_thermal::range& r)
{
  s << "(" << r.min << ", "<<r.max << ")";
  return s;
}

const eos_thermal_impl& eos_thermal::implementation() const
{
  if (pimpl) return *pimpl;
  throw logic_error("EOS: uninitialized use.");
}

eos_thermal_impl::~eos_thermal_impl() {}

/*
This method should be overridden for any EOS that has a more efficient
way of computing this from temperature directly instead of converting
to internal energy first.
*/
pz_real eos_thermal_impl::press_from_valid_rho_temp_ye(const pz_real rho,
            const pz_real temp, const pz_real ye) const
{
  const pz_real eps = eps_from_valid_rho_temp_ye(rho, temp, ye);
  return press_from_valid_rho_eps_ye(rho, eps, ye);
}

/*
This method should be overridden for any EOS that has a more efficient
way of computing this from temperature directly instead of converting
to internal energy first.
*/
pz_real eos_thermal_impl::csnd_from_valid_rho_temp_ye(const pz_real rho,
            const pz_real temp,
            const pz_real ye) const
{
  const pz_real eps = eps_from_valid_rho_temp_ye(rho, temp, ye);
  return csnd_from_valid_rho_eps_ye(rho, eps, ye);
}

/*
This method should be overridden for any EOS that has a more efficient
way of computing this from temperature directly instead of converting
to internal energy first.
*/
pz_real eos_thermal_impl::entropy_from_valid_rho_temp_ye(const pz_real rho,
            const pz_real temp,
            const pz_real ye) const
{
  const pz_real eps = eps_from_valid_rho_temp_ye(rho, temp, ye);
  return entropy_from_valid_rho_eps_ye(rho, eps, ye);
}


pz_real eos_thermal_impl::press_from_rho_eps_ye(const pz_real rho,
            const pz_real eps, const pz_real ye, status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) return nan();
  return press_from_valid_rho_eps_ye(rho, eps, ye);
}

pz_real eos_thermal_impl::csnd_from_rho_eps_ye(const pz_real rho, const pz_real eps,
            const pz_real ye, status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) return nan();
  return csnd_from_valid_rho_eps_ye(rho, eps, ye);
}

pz_real eos_thermal_impl::temp_from_rho_eps_ye(const pz_real rho, const pz_real eps,
            const pz_real ye, status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) return nan();
  return temp_from_valid_rho_eps_ye(rho, eps, ye);
}

void eos_thermal_impl::press_derivs_from_rho_eps_ye(pz_real& press, pz_real& dpdrho,
            pz_real& dpdeps, const pz_real rho, const pz_real eps, const pz_real ye,
            status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) {
    press  = dpdrho = dpdeps = nan();
    return;
  }
  press_derivs_from_valid_rho_eps_ye(press, dpdrho, dpdeps, rho, eps, ye);
}


pz_real eos_thermal_impl::entropy_from_rho_eps_ye(const pz_real rho, const pz_real eps,
            const pz_real ye, status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) return nan();
  return entropy_from_valid_rho_eps_ye(rho, eps, ye);
}

void eos_thermal_impl::press_csnd_from_rho_eps_ye(pz_real& press, pz_real& csnd,
            const pz_real rho, const pz_real eps, const pz_real ye,
            status& stat) const
{
  check_rho_eps_ye(rho, eps, ye, stat);
  if (stat.failed) {
    press  = csnd = nan();
    return;
  }
  press = press_from_valid_rho_eps_ye(rho, eps, ye);
  csnd  = csnd_from_valid_rho_eps_ye(rho, eps, ye);
}

pz_real eos_thermal_impl::eps_from_rho_temp_ye(const pz_real rho, const pz_real temp, const pz_real ye,
            status& stat) const
{
  check_rho_temp_ye(rho, temp, ye, stat);
  if (stat.failed) return nan();
  return eps_from_valid_rho_temp_ye(rho, temp, ye);
}

pz_real eos_thermal_impl::press_from_rho_temp_ye(const pz_real rho,
            const pz_real temp, const pz_real ye, status& stat) const
{
  check_rho_temp_ye(rho, temp, ye, stat);
  if (stat.failed) return nan();
  return press_from_valid_rho_temp_ye(rho, temp, ye);
}

pz_real eos_thermal_impl::csnd_from_rho_temp_ye(const pz_real rho,
            const pz_real temp, const pz_real ye, status& stat) const
{
  check_rho_temp_ye(rho, temp, ye, stat);
  if (stat.failed) return nan();
  return csnd_from_valid_rho_temp_ye(rho, temp, ye);
}

void eos_thermal_impl::press_csnd_from_rho_temp_ye(pz_real& press, pz_real& csnd,
            const pz_real rho, const pz_real temp, const pz_real ye,
            status& stat) const
{
  check_rho_temp_ye(rho, temp, ye, stat);
  if (stat.failed) {
    press  = csnd = nan();
    return;
  }
  const pz_real eps = eps_from_valid_rho_temp_ye(rho, temp, ye);
  press = press_from_valid_rho_eps_ye(rho, eps, ye);
  csnd  = csnd_from_valid_rho_eps_ye(rho, eps, ye);
}

pz_real eos_thermal_impl::entropy_from_rho_temp_ye(const pz_real rho,
            const pz_real temp, const pz_real ye, status& stat) const
{
  check_rho_temp_ye(rho, temp, ye, stat);
  if (stat.failed) return nan();
  return entropy_from_valid_rho_temp_ye(rho, temp, ye);
}


eos_thermal_impl::range eos_thermal_impl::range_eps(const pz_real rho, const pz_real ye) const
{
  if (is_rho_valid(rho) && is_ye_valid(ye))
    return range_eps_from_valid_rho_ye(rho, ye);
  return range(0, 0);
}

pz_real eos_thermal_impl::nan()
{
  return std::numeric_limits<float>::quiet_NaN();
}

void eos_thermal_impl::check_rho_eps_ye(const pz_real rho, const pz_real eps, const pz_real ye, status& stat) const
{
  stat.failed = !is_rho_eps_ye_valid(rho, eps, ye);
  if (stat.failed) {
    stringstream ss;
    ss << "EOS: parameters out of valid range"
       << "(rho = " << rho << ", eps = " << eps << ", ye = " << ye << ")" << endl
       << "Ranges are rho" << range_rho() << ", eps" << range_eps(rho,ye) << ", Ye" << range_ye();
    stat.err_msg = ss.str();
	}
}

void eos_thermal_impl::check_rho_temp_ye(const pz_real rho, const pz_real temp, const pz_real ye, status& stat) const
{
  stat.failed = !is_rho_temp_ye_valid(rho, temp, ye);
  if (stat.failed) {
    stringstream ss;
    ss << "EOS: parameters out of valid range"
       << "(rho = " << rho << ", temp = " << temp << ", ye = " << ye << ")" << endl
       << "Ranges are rho" << range_rho() << ", temp" << range_temp() << ", Ye" << range_ye();
    stat.err_msg = ss.str();
  }
}


