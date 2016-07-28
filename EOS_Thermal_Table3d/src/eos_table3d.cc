#include "eos_tabulated3d_c_api.h"
#include "eos_table3d.h"
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace whizza;

bool eos_table3d::global_eos_initialized = false;
void eos_table3d::init_global_eos(std::string name)
{
	if (global_eos_initialized) {
    throw logic_error("eos_table3d: multiple initialization");
	}
  global_eos_initialized = (0==tablereader(name.c_str()));
	if (!global_eos_initialized) {
    throw logic_error("eos_table3d: initialization failed");
	}
}

eos_table3d::eos_table3d()
{
  if (!global_eos_initialized) {
    throw logic_error("eos_table3d: global eos not initialized");
  }
  range r;
  ::rho_range(&r.min, &r.max);
  set_range_rho(r);
  ::ye_range(&r.min, &r.max);
  set_range_ye(r);
  ::temp_range(&r.min, &r.max);
  set_range_temp(r);
}

eos_table3d::~eos_table3d() {}

double eos_table3d::press_from_valid_rho_eps_ye(const double rho, const double eps, const double ye) const
{
  int ierr(0);
  double press = tab3d_press(rho, eps, ye, &ierr);
  assert(ierr==0); //ierr != 0 on valid input is considered a bug.
  return press;
}

double eos_table3d::csnd_from_valid_rho_eps_ye(const double rho, const double eps, const double ye) const
{
  int ierr(0);
  double csnd2 = tab3d_csnd2(rho, eps, ye, &ierr);
  assert(ierr == 0); //ierr != 0 on valid input is considered a bug.
  assert(csnd2 >= 0); //Soundspeed^2 should never ever be negative
  return sqrt(csnd2);
}

double eos_table3d::temp_from_valid_rho_eps_ye(const double rho, const double eps, const double ye) const
{
  int ierr(0);
  double temp = tab3d_temp(rho, eps, ye, &ierr);
  assert(ierr==0); //ierr != 0 on valid input is considered a bug.
  return temp;
}

void eos_table3d::press_derivs_from_valid_rho_eps_ye(pz_real& press, pz_real& dpdrho, pz_real& dpdeps,
           const pz_real rho, const pz_real eps, const pz_real ye) const
{
  fprintf(stderr, "eos_table3d::press_derivs_from_valid_rho_eps_ye is not supported anymore!");
  exit(EXIT_FAILURE);
}

double eos_table3d::press_from_valid_rho_temp_ye(const double rho, const double temp, const double ye) const
{
    return tab3d_press_from_temp(rho, temp, ye);
}

double eos_table3d::csnd_from_valid_rho_temp_ye(const double rho, const double temp, const double ye) const
{
  double csnd2 = tab3d_csnd2_from_temp(rho, temp, ye);
  assert(csnd2 >= 0); //Soundspeed^2 should never ever be negative
  return sqrt(csnd2);
}

double eos_table3d::entropy_from_valid_rho_temp_ye(const double rho, const double temp, const double ye) const
{
    return tab3d_entropy_from_temp(rho, temp, ye);
}

double eos_table3d::entropy_from_valid_rho_eps_ye(const double rho, const double eps, const double ye) const
{
  int ierr(0);
  double entropy =  tab3d_entropy_from_eps(rho, eps, ye, &ierr);
  assert(ierr==0); //ierr != 0 on valid input is considered a bug.
  return entropy;
}

double eos_table3d::eps_from_valid_rho_temp_ye(const double rho, const double temp, const double ye) const
{
  return tab3d_eps(rho, temp, ye);
}

eos_thermal_impl::range eos_table3d::range_eps_from_valid_rho_ye(const double rho, const double ye) const
{
  range rg;
  int ierr = ::eps_range(rho, ye, &rg.min, &rg.max);
  assert(ierr==0); //ierr != 0 on valid input is considered a bug.
  return rg;
}

eos_thermal whizza::make_eos_table3d()
{
  return eos_thermal(new eos_table3d());
}

