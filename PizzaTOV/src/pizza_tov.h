#include "tovsol.h"
#include "modes.h"
#include "pizzacactus.h"
#include "pizza_idbase.h"
#include <string>

namespace Pizza {
namespace TOV {

using EOS_Barotropic::eos_1p;
using EOS_Barotropic::eos_cold;
using IDBase::idvars_pc;
using IDBase::idvars;

class ptov : cactus_grid {
  idvars_pc idata;
  vec_u_pc shift;
  scalar_pc lapse, temp, efrac;

  eos_1p eos;
  tovsol star;
  const pz_real pert_amp;
  const int pert_l, pert_m, pert_n;
  const pz_real kick_amp;
  bool set_shift, set_lapse, set_temp, set_efrac;

  pz_real perturbation(vec_u x);
  void init_vars(vec_u x, idvars& id, double& id_lapse, vec_u& id_shift);
  void pert_loc(pmode_cowling &mode, const eos_cold& eosc,
                 vec_u p, idvars& id, double scale) const;
  void kick_loc(vec_u p, idvars& id, double scale) const;

  void optional_hydro(const region::iterator& i, const double rmd);
  public:
  ptov(
    const eos_1p& eos_, pz_real rmd_c_, pz_real gm1_cut_, pz_real pert_amp_,
    int pert_l_, int pert_m_, int pert_n_, pz_real kick_amp_, bool set_shift, bool set_lapse,
    bool set_temp_, bool set_efrac_,
    std::string out_dir_
  );

  void prep(const cGH *cgh);
  void initial_data();
  void perturb();
  void kick();
};

}
}

