#include <cstddef>
#include "gsl/gsl_sf_legendre.h"
#include "pizza_tov.h"
#include <algorithm>
#include <fstream>

using namespace Pizza;
using namespace Pizza::TOV;
using namespace std;


ptov::ptov(const eos_1p& eos_, pz_real rmd_c_, pz_real gm1_cut_, pz_real pert_amp_,
    int pert_l_, int pert_m_, int pert_n_, pz_real kick_amp_, bool set_shift_, bool set_lapse_,
    bool set_temp_, bool set_efrac_, string out_dir_)
: eos(eos_), pert_amp(pert_amp_), pert_l(pert_l_), pert_m(pert_m_), pert_n(pert_n_),
kick_amp(kick_amp_),
set_shift(set_shift_), set_lapse(set_lapse_), set_temp(set_temp_), set_efrac(set_efrac_)
{
  star  = tovsol(rmd_c_, eos, 0.0, 1e-10, gm1_cut_);
  if (root()) {
    string fn = out_dir_+"/star_info";
    ofstream os(fn.c_str());
    os << star.to_str();
    fn  = out_dir_+"/star_profile.dat";
    star.save(fn, 800, 2.0);
  }
}


void ptov::init_vars(vec_u x, idvars& id, double& id_lapse, vec_u& id_shift)
{
  const pz_real rc           = x.norm();
  const tov_vars tv          = star(rc);
  const pz_real h            = tv.g_rr();
  const pz_real fbr2         = 2.0 * h * tv.m_by_r3;

  switch (coord_type()) {
    case cartesian:
      error::unless(coord_sym(0)!=flat && coord_sym(1)!=flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");
      id.glo = ONE;
      for (int i=0;i<3;i++) {
        for (int j=0;j<=i;j++) {
          id.glo(i,j) += x(i) * x(j) * fbr2;
        }
      }
      break;
    case cylindrical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");
      id.glo = ONE;
      id.glo(0,0) += x(0) * x(0) * fbr2;
      id.glo(2,0) += x(2) * x(0) * fbr2;
      id.glo(2,2) += x(2) * x(2) * fbr2;
      id.glo(1,1) = sqr(x(0));
      break;
    case spherical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)==flat,
        "TOV: unsupported topology");
      id.glo      = ZERO;
      id.glo(0,0) = h;
      id.glo(1,1) = id.glo(2,2) = sqr(x(0));
      break;
    default:
      assert(false);
  }

  id.klo      = ZERO;
  id_lapse    = tv.lapse();
  id_shift    = ZERO;

  id.rmd      = tv.rmd;
  id.sed      = eos.sed_from_gm1(eos.gm1_from_rmd(tv.rmd));
  id.press    = tv.p;
  id.vel      = ZERO;
  id.florentz = 1.0;
}

void ptov::optional_hydro(const region::iterator& i, const double rmd)
{
  const double gm1 = eos.gm1_from_rmd(rmd);
  if (set_temp)
    temp(i)  = eos.temp_from_gm1(gm1);
  if (set_efrac)
    efrac(i) = eos.efrac_from_gm1(gm1);
}

void ptov::initial_data()
{
  idvars id;
  double id_lapse;
  vec_u id_shift;
  for (region::iterator i(r_all());i;++i) {
    init_vars(coord(i), id, id_lapse, id_shift);

    idata(i)  <<  id;

    if (set_lapse)
      lapse(i) = id_lapse;

    if (set_shift)
      shift(i) << id_shift;

    optional_hydro(i, id.rmd);
  }
}


void ptov::pert_loc(pmode_cowling &mode, const eos_cold& eosc,
                    vec_u p, idvars& id, double scale) const
{
  vec_u er,eth,ephi;
  pz_real r,th,phi;
  coord_spherical(p, r, th, phi);
  pmode_cowling::pert_local pert = mode(r, th, phi);

  switch (coord_type()) {
    case cartesian:
      error::unless(coord_sym(0)!=flat && coord_sym(1)!=flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");

      er(0)   = sin(th)*cos(phi);
      er(1)   = sin(th)*sin(phi);
      er(2)   = cos(th);

      eth(0)  = cos(th)*cos(phi);
      eth(1)  = cos(th)*sin(phi);
      eth(2)  = -sin(th);

      ephi(0) = -sin(phi);
      ephi(1) = cos(phi);
      ephi(2) = 0;

      id.vel  += (pert.dvr.imag() * er + pert.dvth.imag() * eth
                  + pert.dvphi.imag() * ephi)*scale;

      break;
    case cylindrical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");

      er(0)   = sin(th);
      er(1)   = 0;
      er(2)   = cos(th);

      eth(0)  = cos(th);
      eth(1)  = 0;
      eth(2)  = -sin(th);

      id.vel  += (pert.dvr.imag() * er + pert.dvth.imag() * eth)*scale;

      break;
    case spherical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)==flat,
        "TOV: unsupported topology");
      er(0)   = 1;
      er(1)   = 0;
      er(2)   = 0;

      id.vel  += pert.dvr.imag() * er *scale;

      break;
    default:
      assert(false);
  }
  pz_real hm1    = eosc.hm1_from_rmd(id.rmd);
  hm1         = max(0.0, hm1 + pert.dh.imag() * scale);
  id.rmd      = eosc.rmd_from_hm1(hm1);
  id.sed      = eosc.sed_from_hm1(hm1);
  id.press    = eosc.p_from_hm1(hm1);
  error::incase(id.vsqr() >= 0.99999, "TOV: velocity perturbation too large, v>=c");
  id.florentz = id.florentz_from_vel();
}

void ptov::perturb()
{
  if (pert_amp==0) return;

  eos_cold eosc = eos;
  pmode_cowling mode(star, pert_l, pert_m, pert_n, 1e-5);
  const double vmean = sqrt(2.0*mode.kin_energy()/star.grav_mass());
  const double scale = pert_amp / vmean;
  idvars id;
  for (region::iterator i(r_all());i;++i) {
    idata(i) >>  id;
    pert_loc(mode, eosc, coord(i), id, scale);
    idata(i) << id;
    optional_hydro(i, id.rmd);
  }
}

void ptov::kick_loc(vec_u p, idvars& id, double scale) const
{
  pz_real amp = scale / star.radius();
  switch (coord_type()) {
    case cartesian:
      error::unless(coord_sym(0)!=flat && coord_sym(1)!=flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");
      if (p.norm() < star.radius()) {
        id.vel  += p * amp;
      }
      break;
    case cylindrical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)!=flat,
        "TOV: unsupported topology");
      if (sqrt(p(0)*p(0) + p(1)*p(1)) < star.radius()) {
        id.vel(0)  += p(0) * amp;
        id.vel(2)  += p(2) * amp;
      }
      break;
    case spherical:
      error::unless(coord_sym(0)==full && coord_sym(1)==flat && coord_sym(2)==flat,
        "TOV: unsupported topology");
      if (p(0) < star.radius()) {
        id.vel(0)  += p(0)*amp;
      }
      break;
    default:
      assert(false);
  }
  error::incase(id.vsqr() >= 0.99999, "TOV: radial kick too large, v>=c");
  id.florentz = id.florentz_from_vel();
}

void ptov::kick()
{
  if (kick_amp==0) return;

  idvars id;
  for (region::iterator i(r_all());i;++i) {
    idata(i) >>  id;
    kick_loc(coord(i), id, kick_amp);
    idata(i) << id;
  }
}

