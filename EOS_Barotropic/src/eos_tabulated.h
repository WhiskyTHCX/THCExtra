#ifndef EOSTABULATED_H
#define EOSTABULATED_H

#include "eos_barotropic.h"
#include "eos_polytropic.h"
#include "lookuptab.h"
#include <vector>
#include <string>

namespace EOS_Barotropic {


///Tabulated degenerate EOS.
/**
This uses lookup tables with logarithmic interpolation (lookup_table).
The interpolation is not thermodynamically consistent, so
tables with a large number of points (>1000) should be used.
For densities below the smallest tabulated value,
a matching polytropic EOS (see eos_polytrope) is used.
For notation, see eos_cold.
*/
class eos_tabulated : public eos_1p_api {
  Pizza::NumUtils::lookup_table_loglog gm1_rmd, gm1_p, sed_gm1, p_gm1, rmd_gm1, cs2_gm1, hm1_gm1;
  Pizza::NumUtils::lookup_table_logx temp_gm1, efrac_gm1;
  double efrac0, temp0;
  eos_polytrope poly;

  public:

///Constructor
  eos_tabulated(
    const std::vector<double>& gm1,   ///< \f$ g - 1 \f$
    const std::vector<double>& rmd,   ///< Rest mass densities
    const std::vector<double>& sed,   ///< Specific internal energies
    const std::vector<double>& p,     ///< Pressures
    const std::vector<double>& cs2,   ///< Squared sound speeds
    const std::vector<double>& temp,  ///< Temperature if known, else empty vector
    const std::vector<double>& efrac, ///< Electron fraction if known, else empty vector
    bool isentropic_,                 ///< Whether EOS is isentropic
    std::string descr_="Tabulated"    ///< Optional description
  );

  virtual double gm1_from_valid_rmd(const double rmd) const;
  virtual double gm1_from_valid_p(const double p) const;
  virtual double sed_from_valid_gm1(const double gm1) const;
  virtual double ied_from_valid_gm1(const double gm1) const;
  virtual double p_from_valid_gm1(const double gm1) const;
  virtual double rmd_from_valid_gm1(const double gm1) const;
  virtual double hm1_from_valid_gm1(const double gm1) const;
  virtual double csnd2_from_valid_gm1(const double gm1) const;
  virtual double temp_from_valid_gm1(const double gm1) const;
  virtual double efrac_from_valid_gm1(const double gm1) const;
};


}//namespace EOS_Barotropic

#endif

