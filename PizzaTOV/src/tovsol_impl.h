#include "tovsol.h"
#include "splines.h"
#include "odes.h"

namespace Pizza {
namespace TOV {
using NumUtils::cubic_spline;
using NumUtils::opt_const;
using NumUtils::ode_sys_r;

///Implementation belonging to the user interface \ref tovsol.
class tovsol::impl {
  tov_vars cent, surf;
  double gm1_surf;
  double inf_omega;
  double m_inertia;
  double m_bind;
  cubic_spline spl_nu, spl_lambda, spl_ob, spl_dr_ob, spl_rp;
  //double compute_binding_energy() const;
  double rp_vac(const double r) const;
public:
  const eos_1p eos;
  impl(const double rmd_c_, const eos_1p &eos_, const double inf_omega_,
        const double acc, const double gm1_surf_);
  tov_vars operator()(double r) const;
  void save(std::string fn, int n_steps, double size, units u_out) const;
  std::string to_str(units u_out=units::si(),
                     units u_star=units::geom_meter()) const;
  const tov_vars& surface() const {return surf;}
  const tov_vars& center() const {return cent;}
  double grav_mass() const {return surf.m();}
  double binding_energy() const {return m_bind;}
  double baryonic_mass() const {return grav_mass() + binding_energy();}
  double radius() const {return surf.r;}
  double proper_radius() const {return surf.rphys;}
  double omega_inf() const {return inf_omega;}
  double minertia() const {return m_inertia;}
  double rmd_c() const {return cent.rmd;}

};

///TOV equations in a form more suitable for numeric solution.
class tov_ode : public ode_sys_r {
  const eos_1p& eos;
  double dnu_surf, gm1_c, rmd_c, est_rs;
  public:
  enum {DNU=0, LAMBDA=1, OB=2, DR_OB=3, BE=4, RP=5, SIZE=6};
  tov_ode(double rmd_c_, const eos_1p& eos_, double gm1_surf_=0.0);
  double gm1_from_dnu(const double delta_nu) const;
  static double m_by_r3(const double r, const double lambda, const double ed);
  static double dr_wb_by_r(const double r, const double dr_wb, const double rmdh);
  virtual void deriv(const double x, const double y[], double dy[]) const;
  virtual bool stop(const double x, const double y[]) const;
  virtual void init(double& x0, double y0[], double& dx0) const;
  void err_weights(double w_abs[], double& w_rel, double& w_xabs, double& w_xrel) const;
};

}
}

