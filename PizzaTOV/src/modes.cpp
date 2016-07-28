#include "modes_impl.h"
#include "frobenius.h"
#include "spharmonics.h"
#include "unitconv.h"

#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

using namespace Pizza;
using namespace TOV;
using namespace std;
using boost::format;
using NumUtils::spherical_harm;
using NumUtils::vector_harm_radial;
using NumUtils::vector_harm_electric;
using NumUtils::opt_const;

pruessp::pruessp(const tovsol& star_, const int l_) :
  pruess_problem(0.0, star_.radius()),
  star(star_), eos(star.eos()), l(l_)
{
  get_left_coeff();
  get_right_coeff();
}

void pruessp::coefficients(const double x,
        double& p, double& kc, double& kv) const
{
  const tov_vars u  = star(x);
  const double g_rr   = u.g_rr();
  const double lapse  = u.lapse();
  const double lapse2 = pow(lapse,2);
  const double psi    = lapse2 * u.rmd / sqrt(g_rr);
  const double mu     = -2.0 * g_rr * u.m_by_r3;
  const double dlnpsi = (2.0-1.0/u.cs2) * u.dr_nu() - u.dr_lambda();
  const double rf     = pow(x, 2*(l+1));
  p                   = rf * psi;
  kc                  = -l * ( dlnpsi / x + mu * (l+1) );
  kv                  = g_rr / (lapse2 * u.cs2);
}

void pruessp::get_left_coeff()
{
  const tov_vars s  = star.center();
  const double mu     = -2.0 * s.g_rr() * s.m_by_r3;
  vector<double> a(3,0), bc(3,0), bv(3,0);
  a[0]  = 2*(l+1);
  a[1]  = 0;
  a[2]  = 0; //incorrect, but only 0*a[2] will be used
  bc[0] = 0;
  bc[1] = 0;
  bc[2] = l*(2.0-1.0/s.cs2) * s.d2r_nu() + l*l*mu;
  bv[0] = 0;
  bv[1] = 0;
  bv[2] = s.g_rr() /(pow(s.lapse(),2) * s.cs2);

  left_a  = polynomial_r(a);
  left_bc = polynomial_r(bc);
  left_bv = polynomial_r(bv);
}

void pruessp::get_eos_lim_surf(double& n0, double& dhn) const
{
  const double rmd_small = star(star.radius()*(1.0-1e-4)).rmd;
  const double dh = eos.hm1_from_rmd(rmd_small);
  double n[2];
  for (int i=1; i<=2; i++) {
    const double hm1  = i*dh;
    const double cs2  = eos.csnd2_from_hm1(hm1);
    n[i-1]            = hm1 / (cs2*(hm1+1));
  }
  dhn   = (n[1]-n[0]) / dh;
  n0    = n[0] - dhn*dh;
}

void pruessp::get_right_coeff()
{
  vector<double> a(3,0), bc(3,0), bv(3,0);
  const tov_vars s  = star.surface();
  const double lapse2  = pow(s.lapse(),2);
  const double mu     = -2.0 * s.g_rr() * s.m_by_r3;
  double n,dhn;
  get_eos_lim_surf(n,dhn);
  const double k0   = s.dr_nu() * (dhn + n*0.5*(1.0 - s.d2r_nu() / pow(s.dr_nu(),2)));
  const double k1   = -2.0*s.dr_nu() + k0 + s.dr_lambda();
  const double k3   = s.g_rr() / lapse2;
  const double k2   = k3 / s.dr_nu();
  const double dk2  = k3 * (2.0 - 2.0*s.dr_lambda() / s.dr_nu() + s.d2r_nu() / pow(s.dr_nu(), 2));
  a[0]  = n;
  a[1]  = -2*(l+1) / s.r + k1;
  a[2]  = 0;  //incorrect, but only 0*a[2] will be used by Frobenius method
  bc[0] = 0;
  bc[1] = -l * n / s.r;
  bc[2] = (l / s.r) * (-n/s.r - k0 + 2*s.dr_nu() - s.dr_lambda()) + l*(l+1)*mu;
  bv[0] = 0;
  bv[1] = n * k2;
  bv[2] = k2 * k0 + n*dk2;

  right_a  = polynomial_r(a);
  right_bc = polynomial_r(bc);
  right_bv = polynomial_r(bv);
}


void pruessp::left_id(const double ev, const double x, double& u, double& du) const
{
  frobenius_method fm(left_a, left_bc + ev*left_bv);
  polynomial_r sol  = fm.first_solution(0.0);
  polynomial_r der  = sol.deriv();
  u                 = sol(x);
  du                = der(x);
}

void pruessp::right_id(const double ev, const double x, double& u, double& du) const
{
  const frobenius_method fm(right_a, right_bc + ev*right_bv);
  const polynomial_r sol  = fm.first_solution(0.0);
  const polynomial_r der  = sol.deriv();
  const double s          = bnd_right - x;
  u                       = sol(s);
  du                      = -der(s);
}

double pruessp::est_ev(const int n) const
{
  const double dyn_freq = star.center().cs2 / star.radius();
  return pow(2*M_PI*(n+1)*dyn_freq, 2);
}


pmode_cowling::pmode_cowling(const tovsol& star, int l, int m, int n, double acc)
: i(new impl(star, l, m, n, acc)) {}

const pmode_cowling::impl& pmode_cowling::si() const
{
  if (i) return *i;
  throw runtime_error("pmode_cowling: uninitialized use.");
}

pmode_cowling::pert_coeff
pmode_cowling::operator()(const double r) const
{
  return si()(r);
}

pmode_cowling::pert_local
pmode_cowling::operator()(const double r, const double th, const double phi) const
{
  return si()(r,th,phi);
}

int pmode_cowling::l() const {return si().l();}
int pmode_cowling::m() const {return si().m();}
int pmode_cowling::n() const {return si().n();}
double pmode_cowling::freq() const {return si().freq();}
double pmode_cowling::kin_energy() const {return si().kin_energy();}
void pmode_cowling::save_astext(std::string fn, size_t res) const
{
  si().save_astext(fn, res);
}


void check_mode_param(int l, int m, int n)
{
  if (l < 0)
    throw invalid_argument("Oscillations Mode: negative l requested.");
  if (n < 0)
    throw invalid_argument("Oscillations Mode: negative n requested.");
  if ((l==0) && (n==0))
    throw invalid_argument("Oscillations Mode: n=0 for case l=0 requested.");
  if (abs(m) > l)
    throw invalid_argument("Oscillations Mode: |m|>l requested.");
}


double pmode_cowling::freq(const tovsol& star, int l,int n, double acc)
{
  check_mode_param(l, 0, n);
  pruessp pp(star, l);
  const double ev       = pp.eigenvalue(n, acc, 0, pp.est_ev(n));
  return sqrt(ev)/(2*M_PI);
}



pmode_cowling::impl::impl(const tovsol& star_, int l_, int m_, int n_,
                    const double acc):
star(star_),
eos(star_.eos()),
md_l(l_),
md_m(m_),
md_n(n_)
{
  check_mode_param(md_l, md_m, md_n);

  pruessp pp(star, md_l);
  vector<sturm_vars> ef;

  const double ev       = pp.eigenvalue(ef, md_n, acc, 0, pp.est_ev(md_n));
  const double w        = sqrt(ev);
  md_freq               = w / (2*M_PI);

  vector<double> vr, vsc, vrad1, vrad2, vang;

  BOOST_FOREACH(const sturm_vars& i, ef) {
    tov_vars s = star(i.x);
    vr.push_back(s.r);
    vsc.push_back(eos.h_from_rmd(s.rmd) * i.u);
    vrad1.push_back(i.du * s.lapse() / (w * s.g_rr()));
    vrad2.push_back(md_l * i.u * s.lapse() / (w * s.g_rr()));
    vang.push_back(i.u * s.lapse() / w);
  }

  pspl    = ess_pert_spline(vr, vsc, vrad1, vrad2, vang);
  r_surf  = star.radius();
}

double pmode_cowling::impl::kin_energy() const
{
  const size_t nsampl = 200;
  vector<double> xi, yi;
  for (size_t i=0; i<=nsampl; i++) {
    const double x      = (r_surf*i) / nsampl;
    const pert_coeff pc = (*this)(x);
    const tov_vars tv   = star(x);
    const double y      = 0.5 * tv.rmd *x*x* (pow(pc.radial, 2)
                                + md_l*(md_l+1) * pow(pc.angular,2));
    xi.push_back(x);
    yi.push_back(y);
  }
  cubic_spline ed(xi, yi, 0, 0);
  const double ekin  = ed.integral(0, r_surf);
  return ekin;
}

void pmode_cowling::impl::save_astext(string fn, size_t res) const
{
  units u = units::si()/units::geom_meter();
  const double f_si = md_freq / u.freq();
  ofstream o(fn.c_str());
  res++;
  o << format("# l = %d\n# m = %d\n# n = %d\n# f = %.15e Hz\n") % md_l % md_m % md_n % f_si;
  o << format("# %-24s %-24s %-24s %-24s \n") % " r/m" % "scalar" % "radial/c" % "angular/c";

  format fmt("%24.15e %24.15e %24.15e %24.15e \n");
  for (size_t i=0; i <= res; i++) {
    const double x = (r_surf * i) / res;
    pert_coeff pc = (*this)(x);
    o<<fmt % x % pc.scalar % pc.radial % pc.angular;
  }
}

ess_pert_spline::ess_pert_spline(const std::vector<double>& r_, const std::vector<double>& scalar_,
    const std::vector<double>& radial1_, const std::vector<double>& radial2_,
    const std::vector<double>& angular_) :
scalar(r_, scalar_, opt_const::NONE, 0.0),
radial1(r_, radial1_, opt_const::NONE, 0.0),
radial2(r_, radial2_, opt_const::NONE, 0.0),
angular(r_, angular_, opt_const::NONE, 0.0)
{}

ess_pert ess_pert_spline::operator()(double r) const
{
  return ess_pert(scalar(r), radial1(r), radial2(r), angular(r));
}

pmode_cowling::pert_coeff
pmode_cowling::impl::ess2coeff(const ess_pert& ep, const double r) const
{
  if (md_l==0) {
    return pert_coeff(r, ep.scalar, ep.radial1, 0);
  }
  const double s    = r / r_surf;
  const double s1   = pow(s, md_l-1);
  const double s2   = s * s1;
  return pert_coeff(r,
    ep.scalar * s2,
    ep.radial1 * s2 + ep.radial2 * s1 / r_surf,
    ep.angular * s1 / r_surf
  );

}


pmode_cowling::pert_coeff
pmode_cowling::impl::operator()(const double r) const
{
  const ess_pert p  = pspl(r);
  return ess2coeff(p, r);
}

pmode_cowling::pert_local
pmode_cowling::impl::operator()(const double r, const double th, const double phi) const
{
  const pert_coeff pc = (*this)(r);
  return  pert_local(pc, md_l, md_m, th, phi);
}


pmode_cowling::pert_coeff::pert_coeff(double r_, double scalar_, double radial_, double angular_) :
r(r_), scalar(scalar_), radial(radial_), angular(angular_)
{}


pmode_cowling::pert_local::pert_local(const pert_coeff& pc, int l, int m, double th, double phi)
{
  const spherical_harm ylm(l, m, th, phi);
  const vector_harm_radial rlm(ylm);
  const vector_harm_electric psilm(ylm);
  const complex<double> i(0,1);
  dh    = pc.scalar  * ylm.value;
  dvr   = i * pc.radial  * rlm.sc_r;
  dvth  = i * pc.angular * psilm.sc_th;
  dvphi = i * pc.angular * psilm.sc_phi;
}


