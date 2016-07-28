#include "pizzacactus.h"
#include "pizza_idbase.h"
#include <string>

namespace Pizza {
namespace LoreneID {

class loreneid : cactus_grid {
  IDBase::idvars_pc idata;
  scalar_pc lapse, temp, efrac;
  vec_u_pc shift;

  EOS_Barotropic::eos_1p eos;
  bool set_shift, set_lapse, set_temp, set_efrac, do_modify, old_behaviour;
  double mass_conv;
  double rho_cut;
  double scale_vel;
  double scale_above_density;
  double spinup_fac;

  void complete_vars(IDBase::idvars& id, int& crp_cnt);
  void optional_id(const region::iterator& i, IDBase::idvars& id, const double omega);
  void modify_id(const vec_u& x, IDBase::idvars& id, const vec_u& idshift, const double idlapse,
                 const double omega);
  public:
  loreneid(
    const EOS_Barotropic::eos_1p& eos_, bool set_shift_, bool set_lapse_, bool set_temp_, bool set_efrac_,
    bool do_modify_, bool old_behaviour_, double mass_conv_, double rho_cut_, double scale_vel_,
    double spinup_fac_, double scale_above_density_
  );

  void prep(const cGH *cgh);
  void add_missing(int& crap_count, const double omega);

};

}
}
