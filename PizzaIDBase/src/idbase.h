#ifndef IDBASE_H
#define IDBASE_H

#include "pizzacactus.h"
#include "pizza_eos_barotropic.h"

namespace Pizza {
namespace IDBase {

using EOS_Barotropic::eos_1p;
using EOS_Barotropic::eos_cold;

class pizza_idbase_central
{
  static std::auto_ptr<const pizza_idbase_central> sptr;
  public:
  eos_1p eos;
  pizza_idbase_central(const eos_1p& eos_) : eos(eos_) {}
  static void init(const eos_1p& eos_);
  static const pizza_idbase_central& get();
};


struct idvars {
  mats_l glo;
  mats_l klo;
  pz_real rmd,sed,press;
  vec_u vel;
  pz_real florentz;
  pz_real vsqr() const {return glo.contract(vel);}
  pz_real florentz_from_vel() const {return 1.0/sqrt(1.0-vsqr());}
};

struct idvars_pc_ {
  map_cart m;
  typedef idvars base;

  mats_l_pc glo;
  mats_l_pc klo;
  scalar_pc rmd,sed,press;
  vec_u_pc vel;
  scalar_pc florentz;

  void init(const map_cart &m_, pz_real **glo_, pz_real **klo_, pz_real *rmd_,
              pz_real *sed_, pz_real *press_, pz_real **vel_, pz_real *florentz_) {
    m=m_; glo.init(m,glo_); klo.init(m,klo_); rmd.init(m,rmd_);
    sed.init(m,sed_); press.init(m,press_); vel.init(m,vel_);
    florentz.init(m, florentz_);
  }
  template<class F,class B> void apply(F &f,B &b) {
    f(glo,b.glo); f(klo,b.klo); f(rmd,b.rmd);
    f(sed,b.sed); f(press,b.press); f(vel,b.vel); f(florentz, b.florentz);
  }
};

typedef pc_iface< idvars_pc_ > idvars_pc;

}
}
#endif

