#include "eos_barotropic.h"
#include <string>
#include <stdexcept>

using namespace std;
using namespace EOS_Barotropic;

eos_1p_api::~eos_1p_api() {}

double eos_1p_api::hm1_from_valid_rmd(const double rmd) const
{
  return hm1_from_valid_gm1(gm1_from_valid_rmd(rmd));
}

double eos_1p_api::sed_from_valid_rmd(const double rmd) const
{
  return sed_from_valid_gm1(gm1_from_valid_rmd(rmd));
}

double eos_1p_api::p_from_valid_rmd(const double rmd) const
{
  return p_from_valid_gm1(gm1_from_valid_rmd(rmd));
}

double eos_1p_api::csnd2_from_valid_rmd(const double rmd) const
{
  return csnd2_from_valid_gm1(gm1_from_valid_rmd(rmd));
}

double eos_1p_api::temp_from_valid_gm1(const double gm1) const
{
  throw logic_error("EOS: temperature not implemented");
}

double eos_1p_api::efrac_from_valid_gm1(const double gm1) const
{
  throw logic_error("EOS: electron fraction not implemented");
}


double eos_1p_api::sed_from_rmd(const double rmd) const
{
  check_rmd_valid(rmd);
  return sed_from_valid_rmd(rmd);
}

double eos_1p_api::p_from_rmd(const double rmd) const
{
  check_rmd_valid(rmd);
  return p_from_valid_rmd(rmd);
}

double eos_1p_api::csnd2_from_rmd(const double rmd) const
{
  check_rmd_valid(rmd);
  return csnd2_from_valid_rmd(rmd);
}

double eos_1p_api::gm1_from_rmd(const double rmd) const
{
  check_rmd_valid(rmd);
  return gm1_from_valid_rmd(rmd);
}


double eos_1p_api::gm1_from_p(const double p) const
{
  check_p_valid(p);
  return gm1_from_valid_p(p);
}

double eos_1p_api::sed_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return sed_from_valid_gm1(gm1);
}

double eos_1p_api::ied_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return ied_from_valid_gm1(gm1);
}


double eos_1p_api::p_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return p_from_valid_gm1(gm1);
}

double eos_1p_api::rmd_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return rmd_from_valid_gm1(gm1);
}

double eos_1p_api::hm1_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return hm1_from_valid_gm1(gm1);
}

double eos_1p_api::hm1_from_rmd(const double rmd) const
{
  check_rmd_valid(rmd);
  return hm1_from_valid_rmd(rmd);
}

double eos_1p_api::csnd2_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return csnd2_from_valid_gm1(gm1);
}

double eos_1p_api::temp_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return temp_from_valid_gm1(gm1);
}

double eos_1p_api::efrac_from_gm1(const double gm1) const
{
  check_gm1_valid(gm1);
  return efrac_from_valid_gm1(gm1);
}

void eos_1p_api::check_rmd_valid(double rmd) const
{
  if (!rg_rmd.contains(rmd))
    throw runtime_error("EOS: density out of range");
}

void eos_1p_api::check_gm1_valid(double gm1) const
{
  if (!rg_gm1.contains(gm1))
    throw runtime_error("EOS: g-1 out of range");
}

void eos_1p_api::check_p_valid(double p) const
{
  if (!rg_p.contains(p))
    throw runtime_error("EOS: pressure out of range");
}


void eos_1p_api::set_ranges(const range& rg_rmd_)
{
  if (rg_rmd_.xmin < 0) {
    throw logic_error("EOS: minimum valid density has to be >= 0.");
  }
  if (rg_rmd_.xmax <= rg_rmd_.xmin) {
    throw logic_error("EOS: degenerate valid density range.");
  }
  rg_rmd = rg_rmd_;
  rg_gm1 = range(gm1_from_rmd(rg_rmd.xmin), gm1_from_rmd(rg_rmd.xmax));
  rg_p   = range(p_from_gm1(rg_gm1.xmin), p_from_gm1(rg_gm1.xmax));
}

const eos_1p_api& eos_1p::s() const
{
  if (pimpl) return *pimpl;
  throw logic_error("EOS: uninitialized use.");
}


eos_cold::eos_cold(const eos_1p& that) : eos_1p(that)
{
  if (! is_isentropic())
    throw runtime_error("Cold EOS: trying to use non-isentropic EOS as isentropic");
}



