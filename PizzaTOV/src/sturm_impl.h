#include "sturm.h"
#include "functors.h"
#include <memory>
#include <cstring>

namespace Pizza {
namespace TOV {
using NumUtils::function_r2r;
using NumUtils::functor_r2r;

struct pruess_vars {
  double u,pdu;
  pruess_vars() {}
  pruess_vars(double u_, double pdu_) : u(u_), pdu(pdu_) {}
  void scale(const double s) {u*=s; pdu*=s;}
};

class pruess_step {
  double x,p,kc,kv,h;
  public:
  pruess_step(double x_, double p_, double kc_, double kv_, double h_) : x(x_), p(p_), kc(kc_), kv(kv_), h(h_) {}
  void advance(const double ev, const double dx, pruess_vars& v) const;
  void advance(const double ev, pruess_vars& v) const {advance(ev,h,v);}
  void advance(const double ev, pruess_vars& v, sturm_vars& center) const;
  double x_start() const {return x - h/2.0;}
};

class pruess_core
{
  const pruess_problem& problem;
  const int nroots;
  std::vector <pruess_step> grid_left, grid_right;
  double p_mid;

  void set_geom(const std::vector<double>& mesh);
  void margins_id(const double ev, sturm_vars& sl, sturm_vars& sr, pruess_vars& vl, pruess_vars& vr) const;
  void setup_grid(const std::vector<double>& mesh, int i0, int i1, std::vector<pruess_step>& grid);
  double integrate_theta(const pruess_vars& v0, const double ev,
        const std::vector<pruess_step>& g) const;
  void integrate_ef(const sturm_vars& s0, pruess_vars& v, const double ev,
        const std::vector<pruess_step>& g, std::vector<sturm_vars>& ef) const;
  static void bracket_ev(const function_r2r rf,
        double& x0, double& x1, size_t num_exp);
  public:
  pruess_core(const pruess_problem& problem_, int nroots_, const std::vector<double>& mesh);
  double shoot(const double ev) const;
  double eigenvalue(const double acc, double min_ev, double max_ev) const;
  void eigenfunction(const double ev, std::vector<sturm_vars>& ef) const;
};


class pruess_froot : public functor_r2r
{
  const pruess_core& pm;
  public:
  pruess_froot(const pruess_core& pm_) : pm(pm_) {}
  virtual double operator()(double x) const {return pm.shoot(x);}
};

class unigrid_seq
{
  const pruess_problem& problem;
  const int nroots;
  const size_t max_res;
  size_t res,rres;
  std::auto_ptr<const pruess_core> core;

  void set_res(const size_t newres);
  double get_ev(const double acc_root, const double min_ev, const double max_ev) const;
  public:
  unigrid_seq(const pruess_problem& problem_, int nroots_,  size_t min_res_, size_t max_res_);
  double eigenvalue(const double acc, const double min_ev, const double max_ev);
  void eigenfunction(const double ev, std::vector<sturm_vars>& ef) const;
};

}
}

