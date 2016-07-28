#ifndef TOVSOL_H
#define TOVSOL_H
#include <string>

#include "unitconv.h"
#include "eos_barotropic.h"
#include <boost/shared_ptr.hpp>

/*! \mainpage
The library tovlib is used to find solutions of the TOV equations
for arbitrary barotropic equations of state (EOS). For isentropic
EOSs, tovlib can also be used to compute oscillation frequencies and
eigenfunctions in the Cowling approximation (fixed spacetime).

\par Main components
  - \ref Pizza::TOV::tovsol "tovsol" solves the TOV equations for a given central density and EOS.
  - \ref Pizza::TOV::pmode_cowling "pmode_cowling" computes eigenfunctions and frequencies
    of pressure modes in the Cowling approximation, for a given star.
  - \ref Pizza::TOV::tovseq "tovseq" computes solution sequences.
  - Python language binding in a module called pytov.
\par Utilities
  - \ref Pizza::TOV::frobenius_method "frobenius_method" solves singular second order ODE using Frobenius' method.
  - \ref Pizza::TOV::pruess_problem "pruess_problem" solves singular Sturm-Liouville problems using Pruess method together with Frobenius method.
  - \ref Pizza::units "units" converts units.
  - \ref spharmonics.h spherical harmonics and vector spherical harmonics computed in a way that is regular at the axis.

Equations of state are represented by a generic interface
\ref EOS_Barotropic::eos_1p "eos_1p",
with the specialisation \ref EOS_Barotropic::eos_cold "eos_cold"
for isentropic EOSs.
Implementations are included for polytropic EOS
(\ref EOS_Barotropic::eos_polytrope "eos_polytrope")
and tabulated EOS (\ref EOS_Barotropic::eos_tabulated "eos_tabulated").
Both types can be loaded from file in a \ref EOS_Barotropic::load_eos_1p "unified" way.
Custom EOSs can easily be used without modifying the library by either tabulating them
or implementing a custom class derived from the abstract EOS interface.
*/




namespace Pizza {
namespace TOV {

using EOS_Barotropic::eos_1p;

///Set of variables describing a TOV solution at a given point.
/**
These variables locally decribe non-rotating star,
as well as the frame dragging in the slow rotation approximation.
See \ref tovsol and eos_cold for more definitions.
*/
struct tov_vars {
  double r;       ///< circumferential radius
  double rphys;   ///< physical distance to center
  double nu;      ///< \f$ \nu \f$
  double lambda;  ///< \f$ \lambda \f$
  double m_by_r3; ///< \f$ m / r^3 \f$
  double w;       ///< Frame dragging \f$ \omega \f$
  double wb;      ///< \f$ \bar{\omega} \equiv \Omega - \omega \f$
  double dr_wb_by_r;   ///< \f$ (\partial_r \bar{\omega}) / r \f$
  double p;       ///< Pressure \f$ P \f$
  double rmd;     ///< Rest mass density \f$ \rho \f$
  double ied;     ///< Internal energy density \f$ \rho_I=\rho\epsilon \f$
  double cs2;     ///< Soundspeed squared \f$ c_s^2 \f$

/// \f$ \frac{m}{r} \f$
  double m_by_r() const;
/// \f$ m \f$
  double m() const;
/// \f$ \partial_r \bar{\omega} \f$
  double dr_wb() const;
// /// Specific enthalpy \f$ h \f$
  //double h() const;
  //double p() const;
/// Total energy density \f$ \rho_E = (1+\epsilon) \rho \f$
  double ed() const;
/// \f$ \rho h \f$
  double rmdh() const;
/// Sound speed \f$ c_s \f$
  double csnd() const;
/// Lapse function \f$ \alpha = \sqrt{g_{tt}} = e^\nu \f$
  double lapse() const;
/// \f$ g_{rr} = e^{2\lambda} \f$
  double g_rr() const;
/// \f$ g_{tt} = e^{2\nu} \f$
  double g_tt() const;
/// \f$ g_{\theta\theta} = r^2 \f$
  double g_thth() const;
/// \f$ g_{\phi\phi} = r^2 \sin^2(\theta) \f$
  double g_pp(double theta) const;
/// Redshift to infinity \f$ = 1 - \alpha \f$
  double redshift() const;
/// \f$ v^\phi = \bar{\omega} / \alpha \f$, 3-velocity \f$ \phi \f$-component
  double vphi() const;
/// Compute \f$ \partial_r \nu \f$ according to the TOV-equations
  double dr_nu() const;
/// Compute \f$ \partial^2_r \nu \f$ from TOV-equations
  double d2r_nu() const;
/// Compute \f$ \partial_r \lambda \f$ according to the TOV-equations
  double dr_lambda() const;
/// Compute \f$ \partial_r^2 \bar{\omega} \f$ according to [Hartle]
  double d2r_wb() const;
/// Compute derivative of binding energy inside radius.
  double dr_eb() const;
/// Data as text for printing or saving
  std::string to_str(units u=units::si()) const;
/// Format of to_str() as comment string
  static std::string to_str_fmt();
};


///Class representing a solution of the TOV equations describing a nonrotating compact star.
/**
This solves the equilibrium equations for an ideal nonrotating fluid
with a degenerate (cold) EOS in general relativity.
Further, rotational corrections up to first order are computed.
Note this slow rotation approximation only yields a frame dragging coefficient
\f$ \bar{\omega} \f$,
apart from that the solution is the same as for the nonrotating case.
The line element is
\f[
 ds^2 = -e^{2\nu} dt^2 + e^{2\lambda} dr^2
        + r^2  d\theta^2 + r^2 \sin(\theta)^2 (d\phi^2 - 2\omega\, d\phi dt)
\f]
Other quantities involved are
\f[ e^{-2\lambda} = 1 - 2 \frac{m}{r} \f]
\f[ \bar{\omega} = \Omega - \omega \f]
Where \f$ \Omega \f$ denotes the angular velocity seen by a distant observer.
*/
class tovsol {
  class impl;
  boost::shared_ptr<impl> pimpl;
  const impl& s() const;

public:

/// Default constructor. Any attempt to use resulting object throws an exception.
  tovsol() {}
/// Copy constructor
/**
This only copies the pointer to the actual implementation, so that
objects can be copied efficiently by value while the data is stored in memory
only once. Ownership of the pointer is shared, cleanup automatic.
*/
  tovsol(const tovsol& that) : pimpl(that.pimpl) {}

/// Constructor
  tovsol(
    const double rmd_c_,            ///< Central rest mass density \f$ \rho \f$
    const eos_1p &eos_,             ///< The EOS of the star
    const double inf_omega_ = 0.0,  ///< The (formal) rotation rate \f$ \Omega \f$
                                    ///< used when computing the frame dragging.
    const double acc_=1e-10,        ///< Accuracy of ODE integration steps.
    const double gm1_surf_=0.0      ///< Stop integration at hm1=hm1_surf (UGLY!)
  );

/// Assignment operator
  tovsol& operator=(const tovsol& that) {pimpl = that.pimpl; return *this;}

/// Evaluate solution at given radius
  tov_vars operator()(
    double r ///< Circumferential radius \f$ r \f$
  ) const;

/// Save stellar profile as table in text format
  void save(
    std::string fn,         ///< Filename
    int n_steps=100,        ///< Radial resolution in points per stellar radius
    double size=1.0,        ///< Maximum radial coordinate divided by stellar radius
    units u_out=units::si() ///< Unit system used for for saving
  ) const;

/// Description of the stellar parameters as text
  std::string to_str(
    units u_out=units::si(),  ///< Unit system used for description
    units u_star=units::geom_meter() ///< Internal units of the star
  ) const;

/// Solution at surface
  const tov_vars &surface() const;

/// Solution at center
  const tov_vars &center() const;

/// The EOS of the star
  const eos_1p& eos() const;

/// The gravitational mass \f$M\f$ of the star
  double grav_mass() const;

/// Gravitational binding energy
/**
\returns actually a mass, defined as
\f[ M_0 - M = \int_0^R 4 \pi r^2 \rho \left( \sqrt{g_{rr}} -1 -\epsilon \right) \mathrm{d}r \f]
*/
  double binding_energy() const;

/// Total baryonic mass \f$M_0\f$
/**
\returns
\f[ M_0 = \int_0^R 4\pi r^2 \rho \sqrt{g_{rr}} \mathrm{d}r \f]
*/
  double baryonic_mass() const;

/// The circumferential radius of the star
  double radius() const;

/// Proper radius (physical distance center to surface)
  double proper_radius() const;

/// Rotation rate at infinity used only to compute the frame dragging corrections.
  double omega_inf() const;

/// Moment of inertia
  double minertia() const;

/// Central density
  double rmd_c() const;
};

}
}

#endif
