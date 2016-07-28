#include "spharmonics.h"
#include <stdexcept>
#include <cmath>
#include <cstdlib>

using namespace std;
using namespace Pizza;
using namespace NumUtils;

void legendre_poly::zero_l()
{
  l       = 0;
  value   = 1;
  deriv_x = 0;
  prev    = 0;
}

void legendre_poly::next_l()
{
  const double next = ((2*l+1) * x * value - l * prev) / (l+1);
  deriv_x           = (l+1) * value + x * deriv_x;
  prev              = value;
  value             = next;
  l++;
}

legendre_poly::legendre_poly(int l_, double x_) :
x(x_)
{
  if (l_<0)
    throw invalid_argument("legendre poly: undefined for negative l");

  zero_l();
  while (l<l_) next_l();
}

void aslegendre_poly_sph::first_l(int m_)
{
  const int ma  = abs(m_);
  m   = m_;

  b   = 1.0 / sqrt(4.0*M_PI);
  for (int k=1; k<=ma; k++) {
    b *= sqrt(double(2*k+1)/double(2*k));
  }
  bp  = 0;

  if ((m>0) && (1 == (ma%2)))
    b = -b;

  l   = ma;

  if (m==0) {
    deriv_th  = 0;
    value     = b;
  }
  else {
    const double g  = pow(sin_th, ma-1);
    value           = g * sin_th * b;
    deriv_th        = ma * g * cos_th * b;
  }
}

void aslegendre_poly_sph::next_l()
{
  const int ma = abs(m);
  if (l+1 <= ma) {
    if (l+1 == ma) first_l(m);
    else l++;
    return;
  }

  const int lp1 = l+1;
  if (m==0) {
    const double a  = (sqrt(double(4*lp1*lp1 - 1)) / lp1);
    double bn = a * cos_th * b;
    if (l>0)
      bn -= a * bp * l / sqrt(double(4*l*l -1));
    deriv_th        = sqrt((2*l+3) / double(2*l+1))*(cos_th * deriv_th - (l+1) * sin_th * b);
    bp              = b;
    b               = bn;
    value           = b;
  }
  else {
    const double a0 = lp1*lp1 - ma*ma;
    const double a1 = sqrt((4*lp1*lp1 -1) / a0);
    const double a2 = sqrt((l*l - ma*ma) / double(4*l*l -1));
    const double a3 = sqrt(a0 * (2*l+3) / double(2*l+1));
    const double a4 = pow(sin_th, ma-1);
    const double bn = a1 * (cos_th * b - a2 * bp);

    deriv_th        = a4 * ((l+1) * cos_th * bn - a3 * b);
    bp              = b;
    b               = bn;
    value           = sin_th * a4 * b;
  }
  l++;
}



aslegendre_poly_sph::aslegendre_poly_sph(int l_, int m_, double th_) :
l(l_),
m(m_),
th(th_)
{
  if (l<0) {
    throw invalid_argument("assoc. legendre poly: undefined for negative l");
  }
  cos_th = cos(th);
  sin_th = sin(th);

  if (l<abs(m)) {
    value = deriv_th = b = bp = 0;
    return;
  }

  first_l(m);
  while (l < l_) next_l();
}



spherical_harm::spherical_harm(int l, int m, double th, double phi) :
alp(l, m, th)
{
  exp_imphi       = complex<double>(cos(m*phi),sin(m*phi));
  value           = alp.value * exp_imphi;
}

complex<double> spherical_harm::deriv_th() const
{
  return alp.deriv_th * exp_imphi;
}

complex<double> spherical_harm::tm_by_sin_th() const
{
  if (alp.m==0) return 0;
  return (alp.m * alp.b * pow(alp.sin_th, abs(alp.m)-1)) * exp_imphi;
}

void vector_harm_radial::init(const spherical_harm& ylm)
{
  sc_r    = ylm.value;
  sc_th   = 0;
  sc_phi  = 0;
}

void vector_harm_electric::init(const spherical_harm& ylm)
{
  sc_r    = 0;
  sc_th   = ylm.deriv_th();
  const complex<double> i(0,1);
  sc_phi  = i * ylm.tm_by_sin_th();
}

void vector_harm_magnetic::init(const spherical_harm& ylm)
{
  sc_r    = 0;
  sc_phi  = ylm.deriv_th();
  const complex<double> i(0,1);
  sc_th   = -i * ylm.tm_by_sin_th();
}


void ps_vec_harm_mag::init(const spherical_harm& ylm)
{
  vector_harm_magnetic rw(ylm);
  const int l = ylm.alp.l;
  const double f = 1.0 / sqrt(double(l*(l+1)));
  sc_r    = 0.0;
  sc_th   = f * rw.sc_th;
  sc_phi  = f * rw.sc_phi;
}

void ps_vec_harm_elec::init(const spherical_harm& ylm)
{
  vector_harm_electric rw(ylm);
  const int l = ylm.alp.l;
  const double f = 1.0 / sqrt(double(l*(l+1)));
  sc_r    = 0.0;
  sc_th   = f * rw.sc_th;
  sc_phi  = f * rw.sc_phi;
}

