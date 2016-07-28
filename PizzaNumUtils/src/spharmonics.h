#ifndef SPHARMONICS_H
#define SPHARMONICS_H

#include <complex>

/*! \file
Definition of objects representing spherical and vector spherical harmonics
and Legendre polynomials.
*/

namespace Pizza {
namespace NumUtils {

///Legendre polynomial and its derivative.
struct legendre_poly
{
  int l;
  double value, deriv_x, x, prev;
  legendre_poly(int l, double x);
  void zero_l();
  void next_l();
};

///Associated Legendre polynomial and its derivative.
struct aslegendre_poly_sph
{
  int l,m;
  double value, deriv_th, th, cos_th, sin_th, b, bp;
  aslegendre_poly_sph(int l_, int m_, double th_);
  void first_l(int m_);
  void next_l();
};

///Spherical harmonic and its derivative.
struct spherical_harm
{
  aslegendre_poly_sph alp;
  std::complex<double> value, exp_imphi;
  spherical_harm(int l, int m, double th, double phi);
  std::complex<double> deriv_th() const;
  std::complex<double> tm_by_sin_th() const;
  int l() const {return alp.l;}
  int m() const {return alp.m;}
};

///Radial vector harmonic
struct vector_harm_radial
{
  std::complex<double> sc_r, sc_th, sc_phi;

  void init(const spherical_harm& ylm);
  vector_harm_radial(const spherical_harm& ylm) {init(ylm);}
  vector_harm_radial(int l, int m, double th, double phi) {init(spherical_harm(l,m,th,phi));}
};

///Vector spherical harmonic (Regge-Wheeler) of electric type.
struct vector_harm_electric
{
  std::complex<double> sc_r, sc_th, sc_phi;

  void init(const spherical_harm& ylm);
  vector_harm_electric(const spherical_harm& ylm) {init(ylm);}
  vector_harm_electric(int l, int m, double th, double phi) {init(spherical_harm(l,m,th,phi));}
};

///Vector spherical harmonic (Regge-Wheeler) of magnetic type.
struct vector_harm_magnetic
{
  std::complex<double> sc_r, sc_th, sc_phi;

  void init(const spherical_harm& ylm);
  vector_harm_magnetic(const spherical_harm& ylm) {init(ylm);}
  vector_harm_magnetic(int l, int m, double th, double phi) {init(spherical_harm(l,m,th,phi));}
};

///Vector spherical harmonic (Pure spin) of magnetic type.
struct ps_vec_harm_mag
{
  std::complex<double> sc_r, sc_th, sc_phi;

  void init(const spherical_harm& ylm);
  ps_vec_harm_mag(const spherical_harm& ylm) {init(ylm);}
  ps_vec_harm_mag(int l, int m, double th, double phi) {init(spherical_harm(l,m,th,phi));}
};

///Vector spherical harmonic (Pure spin) of electric type.
struct ps_vec_harm_elec
{
  std::complex<double> sc_r, sc_th, sc_phi;

  void init(const spherical_harm& ylm);
  ps_vec_harm_elec(const spherical_harm& ylm) {init(ylm);}
  ps_vec_harm_elec(int l, int m, double th, double phi) {init(spherical_harm(l,m,th,phi));}
};


}
}

#endif

