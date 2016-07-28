#ifndef EOSPOLYTROPIC_H
#define EOSPOLYTROPIC_H

#include "eos_barotropic.h"
#include <string>

namespace EOS_Barotropic {

///Polytropic EOS
/**
The polytropic EOS is written (using units with \f$ c=1 \f$) as
\f[ P = \rho_p \left(\frac{\rho}{\rho_p}\right)^\Gamma ,\qquad
\Gamma= 1+\frac{1}{n}  \f]
Here we use a <em>polytropic density scale</em> \f$ \rho_p \f$ to specify
the EOS instead of the usual form \f$ P = K \rho^\Gamma \f$
because it has simpler units than
\f$ K = \rho_p^{-1/n} \f$.

See eos_cold for notation used and eos_cold_api for a description of the member functions.
*/
class eos_polytrope : public eos_1p_api {
  double n;       ///< Polytropic index \f$ n \f$
  double rmd_p;   ///< Polytropic density scale \f$ \rho_p \f$
  double np1;     ///< \f$ n+1 \f$
  double gamma;   ///< Polytropic exponent \f$ \Gamma \f$
  double invn;    ///< \f$ \frac{1}{n} \f$

  void init(
    double n_,                          ///<Adiabatic index \f$ n \f$
    double rmd_p_,                      ///<Density scale \f$ \rho_p \f$
    double rmd_max_                     ///<Max valid density
  );

  public:

///Constructor
  eos_polytrope(
    double n_,                          ///<Adiabatic index \f$ n \f$
    double rmd_p_,                      ///<Density scale \f$ \rho_p \f$
    double rmd_max_,                    ///<Max valid density
    std::string descr_ = "Polytropic"   ///<Optional short description
  );

///Alternative constructor
  eos_polytrope(
    double rmd_m, ///<Rest mass density \f$ \rho \f$ at matching point
    double sed_m, ///<Specific internal energy \f$ \epsilon \f$ at matching point
    double p_m,   ///<Pressure \f$ P \f$ at matching point
    double rmd_max_,                    ///<Max valid density
    std::string descr_ = "Polytropic" ///<Optional short description
  );


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

