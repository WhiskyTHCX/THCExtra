#ifndef STURM_H
#define STURM_H

#include <vector>
#include <cstring>

namespace Pizza {
namespace TOV {

struct sturm_vars
{
  double x,u,du;
  sturm_vars() {}
  sturm_vars(double x_, double u_, double du_) : x(x_), u(u_), du(du_) {}
  sturm_vars& operator*=(const double s) {u*=s; du*=s; return *this;}
};

///Abstract interface describing a Sturm-Liouville problem and Pruess method to solve it.
class pruess_problem
{
  public:
  const double bnd_left;
  const double bnd_right;
  pruess_problem(double bnd_left_, double bnd_right_) : bnd_left(bnd_left_), bnd_right(bnd_right_) {}
  virtual ~pruess_problem() {}
  virtual void coefficients(const double x, double& p, double& kc,
        double& kv) const=0;
  virtual void left_id(const double ev, const double x, double& u, double& du) const=0;
  virtual void right_id(const double ev, const double x, double& u, double& du) const=0;
  sturm_vars left_id(const double ev, const double x) const;
  sturm_vars right_id(const double ev, const double x) const;
  double eigenvalue(const size_t nroot, const double acc, const double min_ev, const double max_ev) const;
  double eigenvalue(std::vector<sturm_vars>& ef, const size_t nroot, const double acc, const double min_ev,
        const double max_ev) const;
};

}
}

#endif




