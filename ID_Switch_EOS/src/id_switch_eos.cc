#include "id_switch_eos.h"
#include <stdexcept>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace Pizza;
//using namespace Pizza::Base;
//using namespace Pizza::IDBase;
using namespace Pizza::ID_Switch_EOS;


idswitcheos::idswitcheos(const whizza::eos_thermal& eos_,
       bool sync_eps_temp_, bool temp_from_eps_, bool limit_efrac_,
       bool set_entropy_)
: eos(eos_), sync_eps_temp(sync_eps_temp_),
  temp_from_eps(temp_from_eps_), limit_efrac(limit_efrac_),
  set_entropy(set_entropy_)
{}

void check_fail(whizza::eos_thermal::status& stat)
{
  if (stat.failed) {
    stringstream ss;
    ss << "ID_Switch_EOS: initial data incompatible with evolution EOS ("
       << stat.err_msg << ").";
    throw runtime_error(ss.str());
  }
}

void idswitcheos::adapt()
{
  whizza::eos_thermal::status stat;
  double rmd_min = eos.range_rho().min;

  for (region::iterator i(r_all());i;++i) {
    if (gf_rmd(i) <= rmd_min) {
      gf_rmd(i)   = 0.0;
      gf_sed(i)   = 0.0;
      gf_press(i) = 0.0;
      if (sync_eps_temp) {
        gf_temp(i) = 0;
      }
      if (set_entropy) {
        gf_entropy(i) = 0;
      }
    }
    else {
      if (limit_efrac && !eos.is_ye_valid(gf_efrac(i))) {
        whizza::eos_thermal::range rgye = eos.range_ye();
        gf_efrac(i) = min(max(gf_efrac(i), rgye.min), rgye.max);
      }
      double sed_min = eos.range_eps(gf_rmd(i), gf_efrac(i)).min;
      if (sync_eps_temp) {
        if (temp_from_eps) {
          gf_sed(i)   = max(gf_sed(i), sed_min);
          gf_temp(i)  = eos.temp_from_rho_eps_ye(gf_rmd(i), gf_sed(i),
                                                 gf_efrac(i), stat);
          check_fail(stat);
        }
        else {
          double temp_min = eos.range_temp().min;
          gf_temp(i)  = max(gf_temp(i), temp_min);
          gf_sed(i)   = eos.eps_from_rho_temp_ye(gf_rmd(i), gf_temp(i),
                                                 gf_efrac(i), stat);
          check_fail(stat);
          gf_sed(i)   = max(gf_sed(i), sed_min); //Avoid corner case trouble
        }
      }
      else {
        gf_sed(i)   = max(gf_sed(i), sed_min);
      }
      if (set_entropy) {
        gf_entropy(i) = eos.entropy_from_rho_eps_ye(gf_rmd(i), gf_sed(i),
                                                gf_efrac(i), stat);
        check_fail(stat);
      }
      gf_press(i) = eos.press_from_rho_eps_ye(gf_rmd(i), gf_sed(i),
                                              gf_efrac(i), stat);
      check_fail(stat);
    }
  }
}


