#ifndef EOSPIECEWISEPOLY_H
#define EOSPIECEWISEPOLY_H

#include "eos_barotropic.h"
#include <string>
#include <vector>

namespace EOS_Barotropic {


///Class representing a segment of the piecewise polytropic EOS
class eos_poly_piece {
  public:
  double rmd0;
  double dsed;
  double gamma;
  double rmd_p;
  double n;
  double np1;
  double invn;
  double gm10;
  double p0;

  eos_poly_piece() {};
  eos_poly_piece(double rmd0_, double sed0_,
                 double gamma_, double rmd_p_);


  double sed_from_rmd(const double rmd) const;
  double p_from_rmd(const double rmd) const;
  double csnd2_from_rmd(const double rmd) const;

  double gm1_from_rmd(const double rmd) const;
  double gm1_from_p(const double p) const;
  double sed_from_gm1(const double gm1) const;
  double ied_from_gm1(const double gm1) const;
  double p_from_gm1(const double gm1) const;
  double rmd_from_gm1(const double gm1) const;
  double hm1_from_gm1(const double gm1) const;
  double csnd2_from_gm1(const double gm1) const;

};



///Piecewise Polytropic EOS
class eos_piecewise_poly : public eos_1p_api {
  ///The polytropic segments
  std::vector<eos_poly_piece> segments;

  ///Find the segment responsible for a given mass density
  const eos_poly_piece& segment_for_rmd(double rmd) const;
  ///Find the segment responsible for a given mass pressure
  const eos_poly_piece& segment_for_p(double p) const;
  ///Find the segment responsible for a given \f$ g - 1\f$
  const eos_poly_piece& segment_for_gm1(double gm1) const;

  public:

  ///Constructor
  eos_piecewise_poly(
    double rmdp0,                           ///<First segment polytropic density scale
    const std::vector<double>& segm_bound,  ///<Densities of segment boundaries
    const std::vector<double>& segm_gamma,  ///<Segment gammas
    double rmd_max_,                        ///<EOS max valid density
    std::string descr_ = "PiecewisePoly"    ///<Optional short description
  );


  //virtual double sed_from_valid_rmd(const double rmd) const;
  //virtual double p_from_valid_rmd(const double rmd) const;
  //virtual double csnd2_from_valid_rmd(const double rmd) const;

  virtual double gm1_from_valid_rmd(const double rmd) const;
  virtual double gm1_from_valid_p(const double p) const;
  virtual double sed_from_valid_gm1(const double gm1) const;
  virtual double ied_from_valid_gm1(const double gm1) const;
  virtual double p_from_valid_gm1(const double gm1) const;
  virtual double rmd_from_valid_gm1(const double gm1) const;
  virtual double hm1_from_valid_gm1(const double gm1) const;
  virtual double csnd2_from_valid_gm1(const double gm1) const;
};

}//namespace EOS_Barotropic


#endif

