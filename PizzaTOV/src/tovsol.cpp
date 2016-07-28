#include "tovsol_impl.h"

#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_math.h>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace Pizza;
using namespace TOV;

/**
For details on the method see \ref tovsol::impl::impl.
*/
tovsol::tovsol(const double rmd_c_, const eos_1p &eos_, const double inf_omega_,
      const double acc_, const double gm1_surf_) :
pimpl(new impl(rmd_c_, eos_, inf_omega_, acc_, gm1_surf_))
{}

const tovsol::impl& tovsol::s() const
{
  if (pimpl) return *pimpl;
  throw logic_error("TOV Solution: uninitialized use.");
}

tov_vars tovsol::operator()(double r) const {return s()(r);}

void tovsol::save(string fn, int n_steps, double size, units u_out) const
{
  s().save(fn, n_steps, size, u_out);
}

std::string tovsol::to_str(units u_out, units u_star) const
{
  return s().to_str(u_out, u_star);
}

const tov_vars& tovsol::surface() const {return s().surface();}
const tov_vars& tovsol::center() const {return s().center();}
const eos_1p& tovsol::eos() const {return s().eos;}
double tovsol::grav_mass() const {return s().grav_mass();}
double tovsol::binding_energy() const {return s().binding_energy();}
double tovsol::baryonic_mass() const {return s().baryonic_mass();}
double tovsol::radius() const {return s().radius();}
double tovsol::proper_radius() const {return s().proper_radius();}
double tovsol::omega_inf() const {return s().omega_inf();}
double tovsol::minertia() const {return s().minertia();}
double tovsol::rmd_c() const {return s().rmd_c();}



double tov_vars::m_by_r() const {
  return m_by_r3*r*r;
}
double tov_vars::m() const {
  return m_by_r3*gsl_pow_3(r);
}
double tov_vars::dr_wb() const {
  return r*dr_wb_by_r;
}

double tov_vars::ed() const {
  return rmd + ied;
}
double tov_vars::rmdh() const {
  return ed() + p;
}
double tov_vars::csnd() const {
  return sqrt(cs2);
}
double tov_vars::lapse() const {
  return exp(nu);
}
double tov_vars::g_rr() const {
  return exp(2.0*lambda);
}
double tov_vars::g_tt() const {
  return -exp(2.0*nu);
}
double tov_vars::g_thth() const {
  return r*r;
}
double tov_vars::g_pp(double theta) const {
  return pow(r*sin(theta),2);
}
double tov_vars::redshift() const {
  return -gsl_expm1(nu);
}
double tov_vars::vphi() const {
  return wb  / lapse();
}
double tov_vars::dr_nu() const {
  return g_rr()*r*(4.0*M_PI * p  + m_by_r3);
}
double tov_vars::d2r_nu() const {
  return g_rr()*((1.0+2*r*dr_lambda())*(4*M_PI*p+m_by_r3)
          - 4*M_PI*r*rmdh()*dr_nu() - 3.0*m_by_r3 + 4*M_PI*ed());
}
double tov_vars::dr_lambda() const {
  return g_rr()*r*(4.0*M_PI * ed() - m_by_r3 );
}
double tov_vars::d2r_wb() const {
  return 4.0*M_PI*g_rr()*rmdh()*(r*r*dr_wb_by_r + 4.0*wb) -4.0*dr_wb_by_r;
}
double tov_vars::dr_eb() const
{
  return 4.0 * M_PI * pow(r, 2) * (rmd * gsl_expm1(lambda) - ied);
}



tov_ode::tov_ode(double rmd_c_, const eos_1p& eos_, double gm1_surf_) :
ode_sys_r(SIZE),
eos(eos_), rmd_c(rmd_c_)
{
  gm1_c     = eos.gm1_from_rmd(rmd_c_);
  dnu_surf  = gsl_log1p(gm1_c) - gsl_log1p(gm1_surf_);
  est_rs    = sqrt(1.0/(3.0*M_PI*rmd_c_));
}

double tov_ode::gm1_from_dnu(const double delta_nu) const
{
  double gm1 = gm1_c + (1.0 + gm1_c) * gsl_expm1(-delta_nu);
  return max(gm1, 0.0);
}

double tov_ode::m_by_r3(const double r, const double lambda, const double ed)
{
  if (r>0)
    return -0.5 * gsl_expm1(-2.0 * lambda) / (r*r);
  return (4.0/3.0) * M_PI * ed;
}

double tov_ode::dr_wb_by_r(const double r, const double dr_wb, const double rmdhwb)
{
  if (r>0)
    return dr_wb / r;
  return (16.0/5.0) * M_PI * rmdhwb;
}

void tov_ode::deriv(const double r, const double y[], double dy[]) const
{
  tov_vars v;
  v.r           = r;
  v.lambda      = y[LAMBDA];
  v.wb          = y[OB];
  const double gm1 = gm1_from_dnu( y[DNU] );
  v.p           = eos.p_from_gm1(gm1);
  v.rmd         = eos.rmd_from_gm1(gm1);
  v.ied         = eos.ied_from_gm1(gm1);
  v.m_by_r3     = m_by_r3(r, v.lambda, v.ed());
  v.dr_wb_by_r  = dr_wb_by_r(r, y[DR_OB], v.rmdh()*y[OB]);


//v.nu, v.w, v.eb, v.rp not needed here

  dy[DNU]       = v.dr_nu();
  dy[LAMBDA]    = v.dr_lambda();
  dy[OB]        = r * v.dr_wb_by_r;
  dy[DR_OB]     = v.d2r_wb();
  dy[BE]        = v.dr_eb();
  dy[RP]        = sqrt(v.g_rr());
}

void tov_ode::init(double& x0, double y0[], double& dx0) const
{
  x0         = 0;
  dx0        = est_rs / 5000.0;
  y0[DNU]    = 0.0;
  y0[LAMBDA] = 0.0;
  y0[DR_OB]  = 0.0;
  y0[OB]     = 1.0;
  y0[BE]     = 0.0;
  y0[RP]     = 0.0;
}

bool tov_ode::stop(const double x, const double y[]) const
{
  return (y[DNU] > dnu_surf);
}

void tov_ode::err_weights(double w_abs[], double& w_rel, double& w_xabs, double& w_xrel) const
{
  w_abs[DNU]    = 1.0;
  w_abs[LAMBDA] = 1.0;
  w_abs[DR_OB]  = 1.0;
  w_abs[OB]     = 1.0;
  w_abs[BE]     = (4.0/3.0)*M_PI * pow(est_rs,3) * rmd_c;
  w_abs[RP]     = est_rs;
  w_rel         = 0.0;
  w_xabs        = 0.0;
  w_xrel        = 1.0;
}


/**
Solves the TOV equations in the form described by \ref tov_ode,
using the ODE solver \ref tov_math::ode_sys_r.
Simultaneously,
the ODE for the frame dragging coefficient for a slowly rotating star (Hartle)
is solved.
The equtions are formulated in a way that does not degenerate at the surface,
and the radius is determined accurately.
*/
tovsol::impl::impl(const double rmd_c_, const eos_1p &eos_, const double inf_omega_,
                    const double acc, const double gm1_surf_)
: inf_omega(inf_omega_), eos(eos_)
{
  double gm1_min  = eos.range_gm1().xmin * (1.0 + 1e-10);
  gm1_surf        = max(gm1_surf_, gm1_min);
  //cout<< gm1_min<<",  "<<gm1_surf<<",  "<<gm1_surf_<<endl;;
  if (rmd_c_ <= 0)
    throw invalid_argument("TOV Solution: negative central density requested.");
  if (acc <= 0)
    throw invalid_argument("TOV Solution: requested accuracy <= 0");
  if (gm1_surf < 0)
    throw invalid_argument("TOV Solution: specified surface by g = g_surf < 1");
  if (eos.rmd_from_gm1(gm1_surf) >= rmd_c_)
    throw invalid_argument("TOV Solution: specified surface above central density.");


  tov_ode ode(rmd_c_, eos, gm1_surf);

  vector<double> res_r;
  vector<vector<double> > res;

  ode.integrate(acc, res_r, res, 1000000);

  surf.r            = res_r.back();
  surf.rphys        = res[tov_ode::RP].back();
  m_inertia         = (1.0/6.0) * gsl_pow_4(surf.r) /
    (surf.r / 3.0 + res[tov_ode::OB].back() / res[tov_ode::DR_OB].back() );

  surf.lambda       = res[tov_ode::LAMBDA].back();
  surf.m_by_r3      = tov_ode::m_by_r3(surf.r, surf.lambda, 0);
  surf.nu           = -surf.lambda;
  surf.w            = inf_omega * 2.0 * m_inertia / gsl_pow_3(surf.r);
  surf.wb           = inf_omega - surf.w;
  surf.dr_wb_by_r   = 6.0 * m_inertia * inf_omega / gsl_pow_5(surf.r);
  surf.p            = eos.p_from_gm1(gm1_surf);
  surf.ied          = eos.ied_from_gm1(gm1_surf);
  surf.rmd          = eos.rmd_from_gm1(gm1_surf);
  surf.cs2          = eos.csnd2_from_gm1(gm1_surf);

  double ddnu = res[tov_ode::DNU].back() - surf.nu;
  BOOST_FOREACH(double& i, res[tov_ode::DNU]) {
    i -= ddnu;
  }

  double wbs = surf.wb / res[tov_ode::OB].back();
  BOOST_FOREACH(double& i, res[tov_ode::OB]) {
    i *= wbs;
  }
  BOOST_FOREACH(double& i, res[tov_ode::DR_OB]) {
    i *= wbs;
  }

  spl_nu      = cubic_spline(res_r, res[tov_ode::DNU], opt_const::NONE, opt_const::NONE);
  spl_lambda  = cubic_spline(res_r, res[tov_ode::LAMBDA], opt_const::NONE, opt_const::NONE);
  spl_ob      = cubic_spline(res_r, res[tov_ode::OB], opt_const::NONE, opt_const::NONE);
  spl_dr_ob   = cubic_spline(res_r, res[tov_ode::DR_OB], opt_const::NONE, opt_const::NONE);
  spl_rp      = cubic_spline(res_r, res[tov_ode::RP], opt_const::NONE, opt_const::NONE);
  cent        = (*this)(0);
  m_bind      = res[tov_ode::BE].back();//compute_binding_energy();
  if (baryonic_mass()<=0) {
    save("failure_profile", 200, 2.0, units::si());
    cout << to_str();
    throw runtime_error("TOV Solution: FAIL ! obtained negative baryon mass");
  }
}

double tovsol::impl::rp_vac(const double r) const
{
  const double  m   = grav_mass();
  const double emla = sqrt(1.0 - 2.0 * m / r);
  return r * emla + m * log(2.0*(r * (emla + 1.0) - m));
}

tov_vars tovsol::impl::operator()(double r) const
{
  if (r<0)
    throw logic_error("TOV Solution: evaluation at r<0 requested.");

  tov_vars v;

  v.r = r;
  if (r < radius()) {
    v.rphys       = spl_rp(r);
    v.nu          = spl_nu(r);
    v.lambda      = spl_lambda(r);
    const double em1 = gsl_expm1(surf.nu - v.nu);
    const double gm1 = em1 + gm1_surf * (1.0 + em1);
    v.p           = eos.p_from_gm1(gm1);
    v.rmd         = eos.rmd_from_gm1(gm1);
    v.ied         = eos.ied_from_gm1(gm1);
    v.cs2         = eos.csnd2_from_gm1(gm1);
    v.wb          = spl_ob(r);
    v.dr_wb_by_r  = tov_ode::dr_wb_by_r(r, spl_dr_ob(r), v.rmdh() * v.wb);
    v.w           = inf_omega - v.wb;
  } else {
    v.rphys       = surf.rphys + rp_vac(v.r) - rp_vac(surf.r);
    v.nu          = 0.5 * gsl_log1p( -2.0 * grav_mass() / r );
    v.lambda      = -v.nu;
    v.w           = inf_omega * 2.0 * m_inertia / gsl_pow_3(r);
    v.wb          = inf_omega - v.w;
    v.dr_wb_by_r  = 6.0 * inf_omega * m_inertia / gsl_pow_5(r);
    v.p           = 0.0;
    v.rmd         = 0.0;
    v.ied         = 0.0;
    v.cs2         = 0.0;
  }
  v.m_by_r3 = tov_ode::m_by_r3(v.r, v.lambda, v.ed());

  return v;
}


string tov_vars::to_str(units u_out) const
{
  units u(units::geom_meter() / u_out);
  boost::format fmt("%24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e %24.15e");
  fmt % (r * u.length())
      % nu
      % lambda
      % (rmd * u.density())
      % (ied * u.density())
      % (p * u.pressure())
      % lapse()
      % g_rr()
      % (wb * u.freq())
      % (dr_wb_by_r * (u.freq() / u.area()) )
      % (w * u.freq())
      % vphi()
      % (m_by_r3 * u.density())
      % (rphys * u.length())
      % sqrt(cs2);

  return fmt.str();
}

string tov_vars::to_str_fmt()
{
  boost::format fmt("%24s %24s %24s %24s %24s %24s %24s %24s %24s %24s %24s %24s %24s %24s %24s ");
  fmt % "# r"
      % "nu"
      % "lambda"
      % "rmd"
      % "ied"
      % "p"
      % "lapse"
      % "g_rr"
      % "omega_bar"
      % "dr_omega_bar/r"
      % "omega"
      % "v^phi"
      % "m/r^3"
      % "r_phys"
      % "csnd/c";
  return fmt.str();
}


void tovsol::impl::save(string fn, int nsteps, double size,  units u_out) const
{
  ofstream f(fn.c_str());
  f << "# Units: " << u_out << endl
    << tov_vars::to_str_fmt() << endl;
  for (double r=0; r< radius()*size; r+= radius()/nsteps )
    f << (*this)(r).to_str(u_out) << endl;
}


string tovsol::impl::to_str(units u_out, units u_star) const
{
  stringstream o;
  units u(u_star / u_out);
  o.precision(15);
  o << scientific
    << "Star model parameters:" << endl << endl
    << "Mass (grav.)      = " << (grav_mass() * u.mass())       << endl
    << "Mass (baryonic)   = " << (baryonic_mass() * u.mass())   << endl
    << "Radius (circumf.) = " << (radius() * u.length())        << endl
    << "Compactness       = " << (2.0*grav_mass() / radius())   << endl
    << "Moment of inertia = " << (minertia() * u.mom_inertia()) << endl
    << "Omega (infinity)  = " << (omega_inf() * u.freq())       << endl
    << "Central" << endl
    << " Restmass density = " << (rmd_c() * u.density())        << endl
    << " Pressure         = " << (cent.p * u.pressure())        << endl
    << " Soundspeed       = " << (cent.csnd() * u.velocity())   << endl
    << " Omega            = " << (cent.wb *u.freq())            << endl
    << " Lapse            = " << cent.lapse()                   << endl
    << "Redshift at surf. = " << surf.redshift()                << endl
    << "in units (" << u_out << ")" << endl
    << endl
    << "EOS: " << eos.to_str() << endl;
  return o.str();
}





