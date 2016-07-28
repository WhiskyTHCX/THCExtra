#include "tovsol.h"
#include "unitconv.h"
#include "splines.h"
#include "functors.h"
#include <string>

namespace Pizza {
namespace TOV {
/*
///Sequence of TOV solutions
class tovseq {
  const EOS_Barotropic::eos_1p eos;
  double rmd_0,rmd_1;
  Pizza::NumUtils::cubic_spline grav_mass_rmd, radius_rmd, ebind_rmd;
public:
  tovseq(double rmd_0, double rmd_1, const EOS_Barotropic::eos_1p &eos_,
          int n_steps=100, const double gm1_surf_=0.0, const double acc_=1e-10);
  //double rmd_c_from_grav_mass(double m) const;
  void save(std::string fn, Pizza::units u_out=Pizza::units::si()) const;
};
*/

class tov_sequence {
  const EOS_Barotropic::eos_1p eos;
  const double gm1_surf;
  //const double rmd_0,rmd_1;
  const Pizza::NumUtils::interval_closed rmd_range;
  const double acc_star;
  const double acc_root;
  double rmdc_max_mb;  ///< Central density of model with maximum baryon mass.
  double rmdc_max_mg;  ///< Central density of model with maximum grav. mass.
  double max_mb;
  double max_mg;

  double rmdc_from_bmass(double mb, double rmd0, double rmd1) const;
  double rmdc_from_gmass(double mg, double rmd0, double rmd1) const;
  std::string fmt_star(const tovsol& star, Pizza::units uconv) const;
public:
  tov_sequence(double rmd_0_, double rmd_1_, const EOS_Barotropic::eos_1p &eos_,
          double gm1_surf_=0.0, int n_steps_=500,
          double acc_star_=1e-10, double acc_max_=1e-5, double acc_root_=1e-8);
  tovsol get_star(double rhoc) const;
  tovsol operator()(double rhoc) const;
  const Pizza::NumUtils::interval_closed& range_rmdc() const {return rmd_range;}
  double rmdc_max_bmass() const {return rmdc_max_mb;}
  double rmdc_max_gmass() const {return rmdc_max_mg;}
  double max_bmass() const {return max_mb;}
  double max_gmass() const {return max_mg;}
  double bmass_from_rmdc(double rmdc) const;
  double gmass_from_rmdc(double rmdc) const;
  double rmdc_from_bmass_stable(double mb) const;
  double rmdc_from_bmass_unstable(double mb) const;
  double rmdc_from_gmass_stable(double mg) const;
  double rmdc_from_gmass_unstable(double mg) const;
  void save(std::string fn, std::string fn_stable, std::string fn_unstable,
            std::string fn_max_mg, std::string fn_max_mb, int ndivisions,
            Pizza::units uout=Pizza::units::si(),
            Pizza::units useq=Pizza::units::geom_meter()) const;
};

}
}

