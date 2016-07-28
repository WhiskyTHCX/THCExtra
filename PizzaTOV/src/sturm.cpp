#include "sturm_impl.h"
#include "roots.h"
#include <cmath>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <boost/foreach.hpp>
#include <algorithm>
#include <string>
#include <stdexcept>

using namespace std;
using namespace Pizza;
using namespace TOV;

pruess_core::pruess_core(const pruess_problem& problem_, int nroots_, const vector<double>& mesh) :
problem(problem_),
nroots(nroots_)
{
  set_geom(mesh);
}



void pruess_core::set_geom(const vector<double>& mesh)
{
  size_t i_mid        = mesh.size()/2;
  const double x_mid  = mesh[i_mid];
  double kcm, kvm;
  problem.coefficients(x_mid, p_mid, kcm, kvm);
  setup_grid(mesh, 0, i_mid, grid_left);
  setup_grid(mesh, mesh.size()-1, i_mid, grid_right);
}



void pruess_core::setup_grid(const vector<double>& mesh, int i0, int i1, vector<pruess_step>& grid)
{
  grid.clear();
  int dir = (i1>i0) ? 1 : -1;

  for (int i=i0; i != i1; i += dir) {
    const double m  = (mesh[i]+mesh[i+dir])/2.0;
    const double h  = mesh[i+dir] - mesh[i];
    double pm, kcm, kvm;
    problem.coefficients(m, pm, kcm, kvm);
    grid.push_back(pruess_step(m, pm/p_mid, kcm, kvm, h));
  }
}



double pruess_core::integrate_theta(const pruess_vars& v0, const double ev,
    const vector<pruess_step>& g) const
{
  pruess_vars v = v0;
  int nwind     = 0;
  BOOST_FOREACH(const pruess_step& i, g) {
    const double ul = v.u;
    i.advance(ev, v);
    if (v.pdu<0) {
      if ((ul <= 0) && (v.u >= 0) && (v.u > ul)) nwind--;
      else if ((ul >= 0) && (v.u <= 0) && (v.u < ul)) nwind++;
    }
  }
  const double theta  = atan2(v.u, v.pdu) + nwind*2*M_PI;
  return theta;
}

void pruess_core::integrate_ef(const sturm_vars& s0, pruess_vars& v, const double ev,
    const vector<pruess_step>& g, vector<sturm_vars>& ef) const
{
  ef.clear();
  ef.push_back(s0);
  BOOST_FOREACH(const pruess_step& i, g) {
    sturm_vars s;
    i.advance(ev, v, s);
    ef.push_back(s);
  }
}

void pruess_core::margins_id(const double ev, sturm_vars& sl, sturm_vars& sr, pruess_vars& vl, pruess_vars& vr) const
{
  double p, kc, kv;
  sl = problem.left_id(ev, problem.bnd_left);
  sturm_vars slm = problem.left_id(ev, grid_left[0].x_start());
  problem.coefficients(grid_left[0].x_start(), p, kc, kv);
  vl = pruess_vars(slm.u, slm.du * p / p_mid);
  sr = problem.right_id(ev, problem.bnd_right);
  sturm_vars srm = problem.right_id(ev, grid_right[0].x_start());
  problem.coefficients(grid_right[0].x_start(), p, kc, kv);
  vr = pruess_vars(srm.u, srm.du * p / p_mid);
}


double pruess_core::shoot(const double ev) const
{
  sturm_vars sl,sr;
  pruess_vars vl,vr;
  margins_id(ev, sl, sr, vl, vr);
  const double theta_l  = integrate_theta(vl, ev, grid_left);
  const double theta_r  = integrate_theta(vr, ev, grid_right);
  return theta_l - theta_r - M_PI * nroots;
}

void pruess_core::bracket_ev(const function_r2r rf, double& x0, double& x1, size_t num_exp)
{
  const double m = (x0 + x1) / 2.0;

  for (size_t n=0; rf(x0)>0; n++) {
    if (n>=num_exp)
      throw runtime_error("Pruess method: root bracketing failed");
    x0 = m - 2.0*(m-x0);
  }
  for (size_t n=0; rf(x1)<0; n++) {
    if (n>=num_exp)
      throw runtime_error("Pruess method: root bracketing failed");
    x1 = m + 2.0*(x1 - m);
  }
}

double pruess_core::eigenvalue(const double acc, double min_ev, double max_ev) const
{
  function_r2r rf(new pruess_froot(*this));
  bracket_ev(rf, min_ev, max_ev, 20);
  return findroot(rf, min_ev, max_ev, 0, acc);
}

void pruess_core::eigenfunction(const double ev, vector<sturm_vars>& ef) const
{
  vector<sturm_vars> ef_left, ef_right;
  sturm_vars sl,sr;
  pruess_vars vl,vr;

  margins_id(ev, sl, sr, vl, vr);
  integrate_ef(sl, vl, ev, grid_left, ef_left);
  integrate_ef(sr, vr, ev, grid_right, ef_right);

  const double w          = (problem.bnd_right-problem.bnd_left);
  const complex<double> zl(vl.u, w*vl.pdu);
  const complex<double> zr(vr.u, w*vr.pdu);
  const double scale = real(zl/zr);
  reverse(ef_right.begin(), ef_right.end());
  BOOST_FOREACH(sturm_vars& i, ef_right) {
    i *= scale;
  }
  ef.clear();
  ef.insert(ef.end(), ef_left.begin(), ef_left.end());
  ef.insert(ef.end(), ef_right.begin(), ef_right.end());
}


sturm_vars pruess_problem::left_id(const double ev,  const double x) const
{
  double u(0), du(0);
  left_id(ev, x, u, du);
  return sturm_vars(x, u, du);
}

sturm_vars pruess_problem::right_id(const double ev,  const double x) const
{
  double u(0), du(0);
  right_id(ev, x, u, du);
  return sturm_vars(x, u, du);
}



double sinc_pi(const double x)
{
  return gsl_sf_sinc(x/M_PI);
}

double sinhc_pi(const double x)
{
  return (gsl_sf_exprel(x)+gsl_sf_exprel(-x)) / 2.0 ;
}

void pruess_step::advance(const double ev, const double dx, pruess_vars& v) const
{
  const double k    = ev*kv - kc;
  const double w    = sqrt(abs(k));
  const double wdx  = w * dx;

  double du         = 0;
  double dpdu       = 0;

  if (k<0) {
    const double sc = sinhc_pi(wdx);
    const double s  = wdx*sc;
    const double a  = s*s / (cosh(wdx)+1.0);
    du              = v.u * a         + v.pdu * sc * dx / p;
    dpdu            = v.u * w * p * s + v.pdu * a;
  }
  else {
    const double sc = sinc_pi(wdx);
    const double s  = wdx*sc;
    const double a  = -2.0 * pow(sin(wdx/2.0),2);
    du              = v.u * a          + v.pdu * sc * dx / p;
    dpdu            = -v.u * p * w * s + v.pdu * a;
  }
  v.u   += du;
  v.pdu +=dpdu;
}

void pruess_step::advance(const double ev, pruess_vars& v, sturm_vars& center) const
{
  pruess_vars vm = v;
  advance(ev,h/2.0,vm);
  center = sturm_vars(x, vm.u, vm.pdu / p);
  advance(ev,h,v);
}






unigrid_seq::unigrid_seq(const pruess_problem& problem_, int nroots_,  size_t min_res_,
                         size_t max_res_) :
  problem(problem_),
  nroots(nroots_),
  max_res(max_res_)
{
  set_res(min_res_);
}




void unigrid_seq::set_res(const size_t newres)
{
  if (newres>max_res)
    throw runtime_error("Pruess method: excessive refinement requested.");

  const double h      = 2.0 / newres;
  const double hm     = 0.1*h;
  const int nsteps    = ceil((log(2.0-hm) - log(hm)) / (2.0*h));
  if (nsteps<2)
    throw runtime_error("Pruess method: too few grid points requested.");
  vector<double> x;
  for (int k=-nsteps; k<=nsteps; k++) {
    const double sk = tanh(k*h);
    const double xk = ((1.0+sk)*problem.bnd_right + (1.0-sk)*problem.bnd_left) / 2.0;
    x.push_back(xk);
  }

  rres = x.size();
  auto_ptr<const pruess_core> newcore(new pruess_core(problem, nroots, x));
  core = newcore;
  res  = newres;
}

double unigrid_seq::get_ev(const double acc_root, const double min_ev, const double max_ev) const
{
  return  core->eigenvalue(acc_root, min_ev, max_ev);
}

double unigrid_seq::eigenvalue(const double acc, const double min_ev, const double max_ev)
{
  double err;
  const double acc_root = acc / 1e2;
  double ev             = get_ev(acc_root, min_ev, max_ev);

  do {
    set_res(2*res);
    const double nev  = get_ev(acc_root, ev*0.9, ev*1.1);
    err               = abs((nev-ev)/nev);
    ev                = nev;
  } while (err > acc);

  return ev;
}

void unigrid_seq::eigenfunction(const double ev, std::vector<sturm_vars>& ef) const
{
  core->eigenfunction(ev,ef);
}

double pruess_problem::eigenvalue(const size_t nroot, const double acc, const double min_ev, const double max_ev) const
{
  unigrid_seq ugs(*this, nroot,  30, 500000);
  return ugs.eigenvalue(acc, min_ev, max_ev);
}


double pruess_problem::eigenvalue(std::vector<sturm_vars>& ef, const size_t nroot, const double acc, const double min_ev, const double max_ev) const
{
  unigrid_seq ugs(*this, nroot,  30, 500000);
  const double ev = ugs.eigenvalue(acc, min_ev, max_ev);
  ugs.eigenfunction(ev, ef);
  return ev;
}


