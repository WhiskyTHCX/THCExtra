#ifndef EOSBAROTROPIC_H
#define EOSBAROTROPIC_H

#include <boost/shared_ptr.hpp>
#include <string>
#include "intervals.h"

namespace EOS_Barotropic {

///Abstract interface for implementations of degenerate (cold) equations of state.
/**
For notation, definitions, and conventions see eos_1p.
*/
class eos_1p_api {
  public:
  typedef Pizza::NumUtils::interval_closed range;

  private:
  bool isentropic;    ///< Whether EOS is isentropic
  bool temp_avail;    ///< Whether EOS has temperature information
  bool efrac_avail;   ///< Whether EOS has electron fraction information
  std::string descr;  ///< Short description text for EOS.
  range rg_rmd;       ///< Density range in which EOS is valid.
  range rg_gm1;       ///< Range of \f$ g-1 \f$ in which EOS is valid.
  range rg_p;         ///< Pressure range in which EOS is valid.

  protected:
  ///Set valid density range. Has to be called by constructor of each implementation.
  void set_ranges(const range& rg_rmd_);
  ///Check if density is in valid range, else throw exception.
  void check_rmd_valid(double rmd) const;
  ///Check if \f$ g-1 \f$ is in valid range, else throw exception.
  void check_gm1_valid(double gm1) const;
  ///Check if pressure is in valid range, else throw exception.
  void check_p_valid(double p) const;

  public:
  eos_1p_api(bool isentropic_, bool temp_avail_,
             bool efrac_avail_, std::string descr_)
  : isentropic(isentropic_), temp_avail(temp_avail_),
    efrac_avail(efrac_avail_), descr(descr_) {}
  virtual ~eos_1p_api();

///Whether EOS is isentropic
  bool is_isentropic() const {return isentropic;}
///Whether EOS can compute temperature
  bool has_temp() const {return temp_avail;}
///Whether EOS can compute electron fraction
  bool has_efrac() const {return efrac_avail;}
///Returns range of validity for density
  const range& range_rmd() const {return rg_rmd;}
///Returns range of validity for \f$ g-1 \f$
  const range& range_gm1() const {return rg_gm1;}
///Returns range of validity for pressure
  const range& range_p() const {return rg_p;}
///Whether given density is in the valid range of the EOS
  bool is_rmd_valid(double rmd) const {return rg_rmd.contains(rmd);}
///Whether given \f$ g-1 \f$ is in the valid range of the EOS
  bool is_gm1_valid(double gm1) const {return rg_gm1.contains(gm1);}
///Whether given pressure is in the valid range of the EOS
  bool is_p_valid(double p) const {return rg_p.contains(p);}

///Compute \f$ \epsilon \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double sed_from_valid_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute Pressure \f$ P \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double p_from_valid_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double csnd2_from_valid_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute \f$ g-1 \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double gm1_from_valid_rmd(
    const double rmd      ///<Rest mass density  \f$ \rho \f$
  ) const =0;

///Compute \f$ g-1 \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double gm1_from_valid_p(
    const double p        ///<Pressure \f$ P \f$
  ) const=0;

///Compute Specific internal energy \f$\epsilon \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double sed_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const=0;

///Compute internal energy density \f$ \rho_I = \rho \epsilon \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double ied_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const=0;

///Compute Pressure \f$ P \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double p_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const =0;

///Compute Rest mass density \f$ \rho \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double rmd_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const=0;

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double hm1_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const=0;

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double hm1_from_valid_rmd(
    const double rmd      ///< \f$ \rho \f$
  ) const;

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double csnd2_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const=0;

///Compute temperature \f$ T \f$ if available, else throw exception.
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double temp_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute electron fraction \f$ Y_e \f$ if available, else throw exception.
/**Assumes input is in the valid range, no checks are performed.*/
  virtual double efrac_from_valid_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute \f$ \epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double sed_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute Pressure \f$ P \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double p_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double csnd2_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const;

///Compute \f$ g-1 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double gm1_from_rmd(
    const double rmd      ///<Rest mass density  \f$ \rho \f$
  ) const;

///Compute \f$ g-1 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double gm1_from_p(
    const double p        ///<Pressure \f$ P \f$
  ) const;

///Compute Specific internal energy \f$\epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double sed_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute internal energy density \f$ \rho_I = \rho \epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double ied_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute Pressure \f$ P \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double p_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute Rest mass density \f$ \rho \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double rmd_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double hm1_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double hm1_from_rmd(
    const double rmd      ///< \f$ \rho \f$
  ) const;

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double csnd2_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute temperature \f$ T \f$ if available, else throw exception.
/**If input is outside valid region of the EOS, throw exception.*/
  double temp_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;

///Compute electron fraction \f$ Y_e \f$ if available, else throw exception.
/**If input is outside valid region of the EOS, throw exception.*/
  double efrac_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const;



///Short description text
  std::string to_str() const {return descr;}
};


///Class representing generic one-parametric EOS
/**
\par
This holds a pointer with shared ownership to the actual implementation.
Therefore it can be copied around by value without wasting resources
(e.g. in case of tabulated EOS) and without the need to manage resources.
\par Notation
It is assumed that the unit sytem is geometric, i.e. \f$ c=G=1 \f$,
but the exact choice is left to the implementation and the code using it.

<table border>
<tr>
<td> \f$ P \f$      </td><td>Pressure              </td><td><tt>p</tt></td>
</tr><tr>
<td> \f$ \rho_E \f$ </td><td>Total energy density  </td><td><tt>ed</tt></td>
</tr><tr>
<td> \f$ n_B \f$    </td><td>Baryon number density </td><td> </td>
</tr><tr>
<td> \f$ \rho = m_B n_B \f$ </td><td>Formal rest mass density </td><td><tt>rmd</tt></td>
</tr><tr>
<td> \f$ m_B \f$    </td><td>Formal baryon mass    </td><td></td>
</tr><tr>
<td> \f$ \epsilon = \rho_E / \rho -1 \f$ </td>
<td>Specific internal energy</td><td><tt>sed</tt></td>
</tr><tr>
<td> \f$ \rho_I = \rho_E - \rho \f$ </td>
<td>Internal energy density</td><td><tt>ied</tt></td>
</tr><tr>
<td> \f$ h = 1+ \epsilon + \frac{P}{\rho} \f$ </td>
<td>Specific enthalpy including restmass</td>
<td><tt>h</tt></td>
</tr><tr>
<td> \f$ h-1 \f$ </td>
<td>Specific enthalpy excluding restmass</td>
<td><tt>hm1</tt></td>
</tr><tr>
<td>\f$ g = \exp\left(\int_{P(\rho=0)}^P \frac{dP'}{\rho_E(P') + P'}\right) \f$</td>
<td>Hydrostatic potential</td><td><tt>g</tt></td>
</tr><tr>
<td>\f$ g - 1 \f$</td>
<td></td><td><tt>gm1</tt></td>
</tr><tr>
<td>\f$ c_s^2 \f$</td><td>Squared adiabatic soundspeed</td><td><tt>csnd2</tt></td>
</tr>
</table>

\par General EOS properties
The EOSs are parametrized by density and alternatively by the quantity
\f$ g-1 \f$ (which is useful to compute hydrostatic equilibrium).
At zero density, \f$ g-1=0 \f$.
For isentropic EOSs, \f$ g=h \f$.
No assumption is made that the EOS is isentropic, one can have e.g. isothermal
as well. Still, the implementation needs to provide that information.
\par
All implementations have to define a validity range. Calls to the EOS
with input outside the valid range throws an exception.
The mass constant \f$ m_B \f$ is not specified. The EOS only deals with
rest mass density, not baryon number density.
Note the specific internal energy and enthalpy are formal quantities
depending on this choice, only \f$ \rho h, \rho_E \f$ have physical meaning.
*/
class eos_1p {
  boost::shared_ptr<eos_1p_api> pimpl;  ///< Pointer to implementation
  const eos_1p_api& s() const;
  public:

  typedef eos_1p_api::range range;

///Default constructor
/**
Creates uninitialized EOS. Any attemt to use resulting object throws an exception.
One can use it after copying another EOS to it.
*/
  eos_1p() {}

///Constructor from pointer to implementation.
/**
The way to use custom EOS implementations.
The object takes ownership of the pointer, do not delete the pointer.
*/
  eos_1p(
    eos_1p_api* eos   ///< Pointer to implementation
  ) : pimpl(eos) {}

///Copy constructor
  eos_1p(const eos_1p& that) : pimpl(that.pimpl) {}
///Assignment operator.
  eos_1p& operator=(const eos_1p& that) {pimpl=that.pimpl; return *this;}

///Whether EOS is isentropic.
  bool is_isentropic() const
  {
    return s().is_isentropic();
  }

///Whether EOS can compute temperature.
  bool has_temp() const
  {
    return s().has_temp();
  }

///Whether EOS can compute electron fraction.
  bool has_efrac() const
  {
    return s().has_efrac();
  }

///Returns range of validity for density
  const range& range_rmd() const {return s().range_rmd();}
///Returns range of validity for \f$ g-1 \f$
  const range& range_gm1() const {return s().range_gm1();}
///Returns range of validity for pressure
  const range& range_p() const {return s().range_p();}
///Whether given density is in the valid range of the EOS
  bool is_rmd_valid(double rmd) const {return s().is_rmd_valid(rmd);}
///Whether given \f$ g-1 \f$ is in the valid range of the EOS
  bool is_gm1_valid(double gm1) const {return s().is_gm1_valid(gm1);}
///Whether given pressure is in the valid range of the EOS
  bool is_p_valid(double p) const {return s().is_p_valid(p);}

///Compute \f$ \epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double sed_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const
  {
    return s().sed_from_rmd(rmd);
  }

///Compute Pressure \f$ P \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double p_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const
  {
    return s().p_from_rmd(rmd);
  }

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double csnd2_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const
  {
    return s().csnd2_from_rmd(rmd);
  }

///Compute \f$ g-1 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double gm1_from_rmd(
    const double rmd                    ///< Rest mass density  \f$ \rho \f$
  ) const
  {
    return s().gm1_from_rmd(rmd);
  }


///Compute \f$ g-1 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double gm1_from_p(
    const double p                      ///<Pressure \f$ P \f$
  ) const
  {
    return s().gm1_from_p(p);
  }

///Compute specific energy \f$\epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double sed_from_gm1(
    const double gm1                    ///< \f$ g-1 \f$
  ) const
  {
    return s().sed_from_gm1(gm1);
  }

///Compute internal energy density \f$\rho_I \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double ied_from_gm1(
    const double gm1                    ///< \f$ g-1 \f$
  ) const
  {
    return s().ied_from_gm1(gm1);
  }

///Compute Pressure \f$ P \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double p_from_gm1(
    const double gm1                    ///< \f$ g-1 \f$
  ) const
  {
    return s().p_from_gm1(gm1);
  }

///Compute Rest mass density \f$ \rho \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double rmd_from_gm1(
    const double gm1                    ///< \f$ g-1 \f$
  ) const
  {
    return s().rmd_from_gm1(gm1);
  }

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double hm1_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const
  {
    return s().hm1_from_gm1(gm1);
  }

///Compute specific enthalpy (including restmass) \f$ h  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double h_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const
  {
    return 1.0 + hm1_from_gm1(gm1);
  }

///Compute specific enthalpy (excluding restmass) \f$ h-1  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double hm1_from_rmd(
    const double rmd      ///< \f$ \rho \f$
  ) const
  {
    return s().hm1_from_rmd(rmd);
  }

///Compute specific enthalpy (including restmass) \f$ h  \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double h_from_rmd(
    const double rmd      ///< \f$ \rho \f$
  ) const
  {
    return 1.0 + hm1_from_rmd(rmd);
  }

///Compute squared adiabatic soundspeed \f$ c_s^2 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double csnd2_from_gm1(
    const double gm1                    ///< \f$ g-1 \f$
  ) const
  {
    return s().csnd2_from_gm1(gm1);
  }

///Compute temperature \f$ T \f$ if available, else throw exception.
/**If input is outside valid region of the EOS, throw exception.*/
  double temp_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const
  {
    return s().temp_from_gm1(gm1);
  }

///Compute electron fraction \f$ Y_e \f$ if available, else throw exception.
/**If input is outside valid region of the EOS, throw exception.*/
  double efrac_from_gm1(
    const double gm1      ///< \f$ g-1 \f$
  ) const
  {
    return s().efrac_from_gm1(gm1);
  }


///Short description text
  std::string to_str() const {return s().to_str();}
};






///Class representing generic isentropic EOS
/**
This specializes the generic one-parametric EOS to the isentropic case.
\par
For isentropic EOS, \f$ g=h \f$, so we parametrize it in terms
of \f$ h-1 \f$ as well (better suited for <tt>PIZZA</tt> evolution code).
\par
All implementations should guarantee that
  - The formal baryon mass constant \f$ m_B \f$ is chosen such that at zero density
    the rest mass density equals the total energy density \f$ \rho_E \f$.
    This implies \f$ \epsilon = 0\f$ at zero rest mass density.
  - At zero density, \f$ h=1\f$.
  - That does not mean the validity range has to start at zero.
*/
class eos_cold : public eos_1p {
  public:
///Default constructor
/**
Creates uninitialized EOS. Any attemt to use resulting object throws an exception.
One can use it after copying another EOS to it.
*/
  eos_cold() {}

///Copy constructor from generic eos_1p. Must be isentropic, else exception is thrown.
  eos_cold(const eos_1p& that);


///Compute Specific energy \f$\epsilon \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double sed_from_hm1(
    const double hm1                    ///< \f$ h-1 \f$
  ) const
  {
    return sed_from_gm1(hm1);
  }

///Compute Pressure \f$ P \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double p_from_hm1(
    const double hm1                    ///< \f$ h-1 \f$
  ) const
  {
    return p_from_gm1(hm1);
  }

///Compute Rest mass density \f$ \rho \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double rmd_from_hm1(
    const double hm1                    ///< \f$ h-1 \f$
  ) const
  {
    return rmd_from_gm1(hm1);
  }

///Compute squared soundspeed \f$ c_s^2 \f$
/**If input is outside valid region of the EOS, throw exception.*/
  double csnd2_from_hm1(
    const double hm1                    ///< \f$ h-1 \f$
  ) const
  {
    return csnd2_from_gm1(hm1);
  }

};

}// namespace EOS_Barotropic


#endif

