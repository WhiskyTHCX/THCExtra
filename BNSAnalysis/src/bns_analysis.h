#include "pizzacactus.h"

namespace Pizza {
namespace BNSAnalysis {


struct bns_ana_in {
  mats_l glo;
  vec_u  shift;
  vec_u  vel;
  vec_l  scons;
  double w_lorentz;
  double lapse;
  double dens;
  double ye;
  double rho;
  double eps;
  double press;
};

struct bns_ana_out {
  double u_t;
  double s_phi;
  double cons_entropy;
  double dens_noatmo;
  double dens_unbnd_geodesic;
  double dens_unbnd_bernoulli;
  double dens_unbnd_garching;
  double entropy_unbnd;
  double ye_unbnd;
};


struct bns_ana_in_pc_ {
  map_cart m;
  typedef bns_ana_in base;

  mats_l_pc glo;
  vec_u_pc shift;
  vec_u_pc vel;
  vec_l_pc scons;
  scalar_pc w_lorentz;
  scalar_pc lapse;
  scalar_pc dens;
  scalar_pc ye;
  scalar_pc rho;
  scalar_pc eps;
  scalar_pc press;


  void init(const map_cart &m_, pz_real **glo_, pz_real **shift_,
            pz_real **vel_, pz_real **scons_, pz_real *w_lorentz_,
            pz_real *lapse_, pz_real *dens_, pz_real *ye_, pz_real *rho_,
            pz_real *eps_, pz_real *press_)
  {
    m=m_; glo.init(m,glo_); shift.init(m,shift_); vel.init(m,vel_);
    scons.init(m,scons_); w_lorentz.init(m,w_lorentz_);
    lapse.init(m,lapse_); dens.init(m,dens_); ye.init(m,ye_);
    rho.init(m,rho_); eps.init(m, eps_); press.init(m, press_);
  }
  template<class F,class B> void apply(F &f,B &b) {
    f(glo,b.glo); f(shift, b.shift); f(vel,b.vel); f(scons, b.scons);
    f(w_lorentz,b.w_lorentz); f(lapse,b.lapse); f(dens,b.dens);
    f(ye, b.ye); f(rho, b.rho); f(eps, b.eps); f(press, b.press);
  }
};
typedef pc_iface< bns_ana_in_pc_ > bns_ana_in_pc;


class bnsanalysis : cactus_grid {
  const bool set_u_t;
  const bool set_s_phi;
  const bool set_cons_entropy;
  const bool set_dens_noatmo;
  const bool set_dens_unbnd;
  const bool set_entropy_unbnd;
  const bool set_ye_unbnd;
  bool need_entropy;
  const double rho_cut;

  bns_ana_in_pc  vars_in;
  scalar_pc p_entropy;
  scalar_pc p_u_t;
  scalar_pc p_s_phi;
  scalar_pc p_cons_entropy;
  scalar_pc p_dens_noatmo;
  scalar_pc p_dens_unbnd_geodesic;
  scalar_pc p_dens_unbnd_bernoulli;
  scalar_pc p_dens_unbnd_garching;
  scalar_pc p_entropy_unbnd;
  scalar_pc p_ye_unbnd;
  void compute_local(const vec_u& pos, const bns_ana_in& in,
                     bns_ana_out& out, const double& entropy);

  public:

  bnsanalysis(bool set_u_t_, bool set_s_phi_, bool set_cons_entropy_,
    bool set_dens_noatmo_, bool set_dens_unbnd_, bool set_entropy_unbnd_,
    bool set_ye_unbnd_, double rho_cut_);
  void prep(const cGH *cgh);
  void compute();
};

}
}
