#ifndef PIZZAODES_H
#define PIZZAODES_H

#include <vector>
#include <cstring>

namespace Pizza {
namespace NumUtils {

///Abstract interface describing a system of first order ODEs, and a method to solve it.
class ode_sys_r
{
  public:
  const size_t size;
  ode_sys_r(size_t size_) : size(size_) {}
  virtual ~ode_sys_r() {}
  ///Derivatives
  virtual void deriv(const double x, const double y[], double dy[]) const=0;
  ///Stop criterion
  virtual bool stop(const double x, const double y[]) const=0;
  ///Set initial data
  virtual void init(double& x0, double y0[], double& dx0) const=0;
  ///Set up step size control
  virtual void err_weights(double w_abs[], double& w_rel, double& w_xabs, double& w_xrel) const=0;
  ///Integrate ODE until stop criterion is fulfilled.
  void integrate(double acc, std::vector<double>& res_x,
    std::vector<std::vector<double> >&res_y, size_t max_steps) const;
};

}
}
#endif

