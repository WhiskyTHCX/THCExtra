#include "pizza_central.h"
#include "bns_analysis.h"
#include <stdexcept>

using namespace Pizza;
using namespace Pizza::Base;
using namespace Pizza::BNSAnalysis;

bnsanalysis::bnsanalysis(bool set_u_t_, bool set_s_phi_,
    bool set_cons_entropy_, bool set_dens_noatmo_, bool set_dens_unbnd_,
    bool set_entropy_unbnd_, bool set_ye_unbnd_, double rho_cut_)
: set_u_t(set_u_t_), set_s_phi(set_s_phi_), set_cons_entropy(set_cons_entropy_),
  set_dens_noatmo(set_dens_noatmo_), set_dens_unbnd(set_dens_unbnd_),
  set_entropy_unbnd(set_entropy_unbnd_), set_ye_unbnd(set_ye_unbnd_),
  rho_cut(rho_cut_)
{
  need_entropy = set_cons_entropy || set_entropy_unbnd;
}


void bnsanalysis::compute_local(const vec_u& pos, const bns_ana_in& in,
                                bns_ana_out& out, const double& entropy)
{
  bool is_atmo      = in.rho <= rho_cut;
  out.dens_noatmo   = is_atmo ? 0 : in.dens;

  double shvel      = in.glo.contract(in.shift, in.vel);
  double pbrho      = in.press / in.rho;
  double h          = 1 + in.eps + pbrho;
  double u_0        = in.w_lorentz * (shvel - in.lapse);

  double c_geo      = -u_0 - 1.0;
  if (c_geo > 0)
    out.dens_unbnd_geodesic = out.dens_noatmo;
  else
    out.dens_unbnd_geodesic = 0.0;

  double c_ber      = -h*u_0 - 1.0;
  if (c_ber > 0)
    out.dens_unbnd_bernoulli = out.dens_noatmo;
  else
    out.dens_unbnd_bernoulli = 0.0;

  double c_gar      = c_ber - (in.lapse / in.w_lorentz) * pbrho;
  if (c_gar > 0)
    out.dens_unbnd_garching = out.dens_noatmo;
  else
    out.dens_unbnd_garching = 0.0;

  out.entropy_unbnd = out.dens_unbnd_geodesic * entropy;
  out.ye_unbnd      = out.dens_unbnd_geodesic * in.ye;
  out.u_t           = u_0;

  out.s_phi         = -pos(1)*in.scons(0) + pos(0)*in.scons(1);
  out.cons_entropy  = out.dens_noatmo * entropy;
}


void bnsanalysis::compute()
{
  bns_ana_in  in;
  bns_ana_out out;

  for (region::iterator i(r_all());i;++i) {
    vars_in(i)   >> in;

    double entropy = need_entropy ? p_entropy(i) : 0.0;

    compute_local(coord(i), in, out, entropy);

    if (set_u_t)           p_u_t(i)           = out.u_t;
    if (set_s_phi)         p_s_phi(i)         = out.s_phi;
    if (set_cons_entropy)  p_cons_entropy(i)  = out.cons_entropy;
    if (set_dens_noatmo)   p_dens_noatmo(i)   = out.dens_noatmo;
    if (set_dens_unbnd) {
      p_dens_unbnd_geodesic(i)    = out.dens_unbnd_geodesic;
      p_dens_unbnd_bernoulli(i)   = out.dens_unbnd_bernoulli;
      p_dens_unbnd_garching(i)   = out.dens_unbnd_garching;
    }
    if (set_entropy_unbnd) p_entropy_unbnd(i) = out.entropy_unbnd;
    if (set_ye_unbnd)      p_ye_unbnd(i)      = out.ye_unbnd;

  }
}


