#include "assert.h"
#include "math.h"
#include "stdlib.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define MAX(A,B) ((A)>(B)?(A):(B))

static
CCTK_REAL  invg4(CCTK_REAL gd[4][4], CCTK_REAL gu[4][4])
{
    // this was done by mathematica by drag and drop
  CCTK_REAL a = gd[0][0];
  CCTK_REAL b = gd[0][1];
  CCTK_REAL c = gd[0][2];
  CCTK_REAL d = gd[0][3];
  CCTK_REAL e = gd[1][1];
  CCTK_REAL f = gd[1][2];
  CCTK_REAL g = gd[1][3];
  CCTK_REAL h = gd[2][2];
  CCTK_REAL i = gd[2][3];
  CCTK_REAL j = gd[3][3];

  CCTK_REAL detg = -(pow(c,2)*pow(g,2)) + a*pow(g,2)*h +
        pow(d,2)*(-pow(f,2) + e*h) + 2*b*c*g*i - 2*a*f*g*i -
        pow(b,2)*pow(i,2) + a*e*pow(i,2) +
        2*d*(c*f*g - b*g*h - c*e*i + b*f*i) + (pow(c,2)*e -
            2*b*c*f + a*pow(f,2) + pow(b,2)*h - a*e*h)*j;
  CCTK_REAL oodetg = 1./detg;

  gu[0][0] = oodetg*(pow(g,2)*h - 2*f*g*i + pow(f,2)*j + e*(pow(i,2) - h*j));
  gu[0][1] = oodetg*(-(d*g*h) + d*f*i + c*g*i - b*pow(i,2) - c*f*j + b*h*j);
  gu[0][2] = oodetg*(d*f*g - c*pow(g,2) - d*e*i + b*g*i + c*e*j - b*f*j);
  gu[0][3] = oodetg*(-(d*pow(f,2)) + c*f*g + d*e*h - b*g*h - c*e*i + b*f*i);
  gu[1][1] = oodetg*(pow(d,2)*h - 2*c*d*i + pow(c,2)*j + a*(pow(i,2) - h*j));
  gu[1][2] = oodetg*(-(pow(d,2)*f) + c*d*g + b*d*i - a*g*i - b*c*j + a*f*j);
  gu[1][3] = oodetg*(c*d*f - pow(c,2)*g - b*d*h + a*g*h + b*c*i - a*f*i);
  gu[2][2] = oodetg*(pow(d,2)*e - 2*b*d*g + pow(b,2)*j + a*(pow(g,2) - e*j));
  gu[2][3] = oodetg*(-(c*d*e) + b*d*f + b*c*g - a*f*g - pow(b,2)*i + a*e*i);
  gu[3][3] = oodetg*(pow(c,2)*e - 2*b*c*f + pow(b,2)*h + a*(pow(f,2) - e*h));

  gu[1][0] = gu[0][1];
  gu[2][0] = gu[0][2];
  gu[3][0] = gu[0][3];
  gu[2][1] = gu[1][2];
  gu[3][1] = gu[1][3];
  gu[3][2] = gu[2][3];

  return detg;
}


static
void    setLAMBDA(CCTK_REAL LAMBDA[4][4], CCTK_REAL LAMBDAi[4][4], CCTK_REAL xix,CCTK_REAL xiy,CCTK_REAL xiz)
{

  // 1e-10 saves the day if no boost is set and you have to divide by 0 inside LAMBDA
  CCTK_REAL xi = sqrt(xix*xix + xiy*xiy + xiz*xiz)+1e-13;
  CCTK_REAL gb = 1./sqrt(1.-xi*xi);

  LAMBDA[0][0] = gb;
  LAMBDA[0][1] = LAMBDA[1][0] = gb*xix;
  LAMBDA[0][2] = LAMBDA[2][0] = gb*xiy;
  LAMBDA[0][3] = LAMBDA[3][0] = gb*xiz;
  LAMBDA[1][1] = (1.+(gb-1.)*(xix*xix)/(xi*xi));
  LAMBDA[2][2] = (1.+(gb-1.)*(xiy*xiy)/(xi*xi));
  LAMBDA[3][3] = (1.+(gb-1.)*(xiz*xiz)/(xi*xi));
  LAMBDA[1][2] = LAMBDA[2][1] = (gb-1.)*xix*xiy/(xi*xi);
  LAMBDA[1][3] = LAMBDA[3][1] = (gb-1.)*xix*xiz/(xi*xi);
  LAMBDA[2][3] = LAMBDA[3][2] = (gb-1.)*xiy*xiz/(xi*xi);

  invg4(LAMBDA,LAMBDAi);
}

static
void    Gamma44(CCTK_REAL g[4][4], CCTK_REAL dg[4][4][4], CCTK_REAL Gamma[4][4][4])
{
  int o,p,q,r;
  CCTK_REAL gi[4][4];

  invg4(g,gi);

  for (o=0; o<=3; o++) {
    for (p=0; p<=3; p++) {
      for (q=0; q<=3; q++) {
        Gamma[o][p][q] = 0.;
        for (r=0; r<=3; r++) {
          Gamma[o][p][q] += 0.5*gi[o][r]* (dg[q][p][r] + dg[p][q][r] - dg[r][p][q]);
        }
      }
    }
  }
}

static
void set_PUNC(int siz,
    CCTK_REAL M,
    CCTK_REAL eps,
    CCTK_REAL px, CCTK_REAL py, CCTK_REAL pz,
    CCTK_REAL mx, CCTK_REAL my, CCTK_REAL mz,
    CCTK_REAL const * xp, CCTK_REAL const * yp, CCTK_REAL const * zp,
    CCTK_REAL * psi4, CCTK_REAL * dpsi4x, CCTK_REAL * dpsi4y, CCTK_REAL * dpsi4z,
    CCTK_REAL * alpha, CCTK_REAL * dalphax, CCTK_REAL * dalphay, CCTK_REAL * dalphaz)
{
  CCTK_REAL r, x,y,z;

  CCTK_REAL LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);

  for(int ijk = 0; ijk < siz; ++ijk) {
    // set transformed coordinates
    x = LAMBDA[1][1]*(xp[ijk]-px) + LAMBDA[1][2]*(yp[ijk]-py) + LAMBDA[1][3]*(zp[ijk]-pz);
    y = LAMBDA[2][1]*(xp[ijk]-px) + LAMBDA[2][2]*(yp[ijk]-py) + LAMBDA[2][3]*(zp[ijk]-pz);
    z = LAMBDA[3][1]*(xp[ijk]-px) + LAMBDA[3][2]*(yp[ijk]-py) + LAMBDA[3][3]*(zp[ijk]-pz);
    r = sqrt(x*x+y*y+z*z);
    r = MAX(r, eps);

    psi4[ijk]    = pow(1.0+0.5*M/r,4);
    dalphax[ijk] = 0; // M/((0.5*M+r)*(0.5*M+r))*x/r;
    dalphay[ijk] = 0; // M/((0.5*M+r)*(0.5*M+r))*y/r;
    dalphaz[ijk] = 0; // M/((0.5*M+r)*(0.5*M+r))*z/r;
    dpsi4x[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*x/r;
    dpsi4y[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*y/r;
    dpsi4z[ijk]  = -2.*M*pow(1.+0.5*M/r,3)/(r*r)*z/r;

    alpha[ijk]  = 1; // (1.0-0.5*M/r)/(1.0+0.5*M/r);
  }
}

static
void boost_spacetime(int siz,
    CCTK_REAL mx, CCTK_REAL my, CCTK_REAL mz,
    CCTK_REAL const * xp, CCTK_REAL const * yp, CCTK_REAL const * zp,
    CCTK_REAL * psi4, CCTK_REAL * dpsi4x, CCTK_REAL * dpsi4y, CCTK_REAL * dpsi4z,
    CCTK_REAL * alpha, CCTK_REAL * dalphax, CCTK_REAL * dalphay, CCTK_REAL * dalphaz,
    CCTK_REAL * betax, CCTK_REAL * betay, CCTK_REAL * betaz,
    CCTK_REAL * gxx, CCTK_REAL * gxy, CCTK_REAL * gxz, CCTK_REAL * gyy, CCTK_REAL * gyz, CCTK_REAL * gzz,
    CCTK_REAL * Kxx, CCTK_REAL * Kxy, CCTK_REAL * Kxz, CCTK_REAL * Kyy, CCTK_REAL * Kyz, CCTK_REAL * Kzz)
{
  int o,p,q,r,s,t;
  CCTK_REAL g[4][4],  u[4],  delg[4][4][4],  Gamma[4][4][4];  // BEFORE BOOST
  CCTK_REAL gp[4][4], up[4], delgp[4][4][4], Gammap[4][4][4]; // AFTER  BOOST (p for primed)
  CCTK_REAL gip[4][4];
  (void)delgp;
#if 0
  CCTK_REAL W,v2;
#endif

  CCTK_REAL * dpsi4[4] = {NULL, dpsi4x, dpsi4y, dpsi4z};
  CCTK_REAL * dalpha[4] = {NULL, dalphax, dalphay, dalphaz};

  CCTK_REAL LAMBDA[4][4],LAMBDAi[4][4];
  setLAMBDA(LAMBDA,LAMBDAi, mx,my,mz);

  for(int ijk = 0; ijk < siz; ++ijk) {
    // do we need to boost?
    // for analytical issues see src/matter/doc/matter.tex
    if (mx!=0. || my!=0. || mz!=0.) {

#if 0
      // set T_munu ... but better to use/transform only u^mu
      // rho p epsl are invariant under transformation
      if (useMATTER) {
        v2 = 0.;
        W  = 1.0/sqrt(1.0 - v2);
        u[0] = W/alpha[ijk];
        u[1] = u[2] = u[3] = 0.;
      }
#endif

      // set g_\mu\nu
      g[0][0] = -(alpha[ijk]*alpha[ijk]);
      for (o=1; o<=3; o++)
        g[0][o] = g[o][0] = 0.;
      for (o=1; o<=3; o++)
      for (p=1; p<=3; p++)
        g[o][p] = psi4[ijk] * (CCTK_REAL)(o==p);

      // time derivative of g_\mu\nu
      for (o=0; o<=3; o++)
      for (p=0; p<=3; p++)
        delg[0][o][p] = 0.;
          // spatial derivative of g_\mu\nu
      for (o=1; o<=3; o++) {
              // lapse
        delg[o][0][0] = -2.*alpha[ijk] * dalpha[o][ijk];
              // shift
        for (p=1; p<=3; p++)
          delg[o][p][0] = delg[o][0][p] = 0.;
              // metric
        for (p=1; p<=3; p++)
        for (q=1; q<=3; q++)
          delg[o][p][q] = dpsi4[o][ijk] * (CCTK_REAL)(q==p);
      }

      // constract the 4 Gamma
      Gamma44(g,delg, Gamma);

      // make coordinate transformation for u^\mu, g_\mu\nu,  Gamma^\sigma_\mu\nu
      for (o=0; o<=3; o++)
      for (p=0; p<=3; p++)
      for (q=0; q<=3; q++) {
        up[o] = 0.;
        gp[o][p] = 0.;
        Gammap[o][p][q] = 0.;
        for (r=0; r<=3; r++) {
          up[o] += LAMBDAi[r][o] * u[r];
          for (s=0; s<=3; s++) {
            gp[o][p] += LAMBDA[o][r] * LAMBDA[p][s] * g[r][s];
            for (t=0; t<=3; t++) {
              Gammap[o][p][q] += LAMBDAi[o][r] * LAMBDA[p][s] * LAMBDA[q][t] * Gamma[r][s][t];
              Gammap[o][p][q] += LAMBDAi[o][r] * 0.; // lucky coincidence ... 0 because Lambda is constant
            }
          }
        }
      }

      // compute inverse gp^{\mu\nu}
      invg4(gp,gip);

      // set transformed g_munu values to grid
      gxx[ijk]    = gp[1][1];
      gxy[ijk]    = gp[1][2];
      gxz[ijk]    = gp[1][3];
      gyy[ijk]    = gp[2][2];
      gyz[ijk]    = gp[2][3];
      gzz[ijk]    = gp[3][3];
      betax[ijk]  = -gip[0][1]/gip[0][0];
      betay[ijk]  = -gip[0][2]/gip[0][0];
      betaz[ijk]  = -gip[0][3]/gip[0][0];
      alpha[ijk]  = sqrt(-1./gip[0][0]);

#if 0
      // set transformed u^mu values to grid
      if (useMATTER) {
        W        = up[0]*alpha[ijk];
        vx[ijk]  = up[1]/W + betax[ijk]/alpha[ijk];
        vy[ijk]  = up[2]/W + betay[ijk]/alpha[ijk];
        vz[ijk]  = up[3]/W + betaz[ijk]/alpha[ijk];
      }
#endif

      // set K_ij
      Kxx[ijk] = -alpha[ijk] * Gammap[0][1][1];
      Kxy[ijk] = -alpha[ijk] * Gammap[0][1][2];
      Kxz[ijk] = -alpha[ijk] * Gammap[0][1][3];
      Kyy[ijk] = -alpha[ijk] * Gammap[0][2][2];
      Kyz[ijk] = -alpha[ijk] * Gammap[0][2][3];
      Kzz[ijk] = -alpha[ijk] * Gammap[0][3][3];
    } else {
      alpha[ijk] = 1.;
      betax[ijk] = 0.;
      betay[ijk] = 0.;
      betaz[ijk] = 0.;

      gxx[ijk] = psi4[ijk];
      gxy[ijk] = 0.;
      gxz[ijk] = 0.;
      gyy[ijk] = psi4[ijk];
      gyz[ijk] = 0.;
      gzz[ijk] = psi4[ijk];

      Kxx[ijk] = 0.;
      Kxy[ijk] = 0.;
      Kxz[ijk] = 0.;
      Kyy[ijk] = 0.;
      Kyz[ijk] = 0.;
      Kzz[ijk] = 0.;
    }
  }
}

void AddPuncture(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  if(verbose) {
      CCTK_INFO("AddPuncture");
  }

  int const siz = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];

  // Puncture spacetime
  CCTK_REAL * psi4_bh    = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dpsi4x_bh  = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dpsi4y_bh  = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dpsi4z_bh  = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * alpha_bh   = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dalphax_bh = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dalphay_bh = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * dalphaz_bh = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * betax_bh   = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * betay_bh   = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * betaz_bh   = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gxx_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gxy_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gxz_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gyy_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gyz_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * gzz_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kxx_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kxy_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kxz_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kyy_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kyz_bh     = malloc(siz*sizeof(CCTK_REAL));
  CCTK_REAL * Kzz_bh     = malloc(siz*sizeof(CCTK_REAL));

  // Check that we have not run out of memory
  assert(psi4_bh);
  assert(dpsi4x_bh);
  assert(dpsi4y_bh);
  assert(dpsi4z_bh);
  assert(alpha_bh);
  assert(dalphax_bh);
  assert(dalphay_bh);
  assert(dalphaz_bh);
  assert(betax_bh);
  assert(betay_bh);
  assert(betaz_bh);
  assert(gxx_bh);
  assert(gxy_bh);
  assert(gxz_bh);
  assert(gyy_bh);
  assert(gyz_bh);
  assert(gzz_bh);
  assert(Kxx_bh);
  assert(Kxy_bh);
  assert(Kxz_bh);
  assert(Kyy_bh);
  assert(Kyz_bh);
  assert(Kzz_bh);

  // Setup a puncture
  set_PUNC(siz,
      puncture_mass,
      puncture_eps,
      puncture_pos_x, puncture_pos_y, puncture_pos_z,
      -puncture_vel_x, -puncture_vel_y, -puncture_vel_z,
      x, y, z,
      psi4_bh, dpsi4x_bh, dpsi4y_bh, dpsi4z_bh,
      alpha_bh, dalphax_bh, dalphay_bh, dalphaz_bh);

  // Boost spacetime
  boost_spacetime(siz,
      -puncture_vel_x, -puncture_vel_y, -puncture_vel_z,
      x, y, z,
      psi4_bh, dpsi4x_bh, dpsi4y_bh, dpsi4z_bh,
      alpha_bh, dalphax_bh, dalphay_bh, dalphaz_bh,
      betax_bh, betay_bh, betaz_bh,
      gxx_bh, gxy_bh, gxz_bh, gyy_bh, gyz_bh, gzz_bh,
      Kxx_bh, Kxy_bh, Kxz_bh, Kyy_bh, Kyz_bh, Kzz_bh);

  // Superimpose puncture spacetime with current spactime
  for(int ijk = 0; ijk < siz; ++ijk) {
    alp[ijk]   += alpha_bh[ijk];
    betax[ijk] += betax_bh[ijk];
    betay[ijk] += betay_bh[ijk];
    betaz[ijk] += betaz_bh[ijk];
    gxx[ijk]   += gxx_bh[ijk];
    gxy[ijk]   += gxy_bh[ijk];
    gxz[ijk]   += gxz_bh[ijk];
    gyy[ijk]   += gyy_bh[ijk];
    gyz[ijk]   += gyz_bh[ijk];
    gzz[ijk]   += gzz_bh[ijk];
    kxx[ijk]   += Kxx_bh[ijk];
    kxy[ijk]   += Kxy_bh[ijk];
    kxz[ijk]   += Kxz_bh[ijk];
    kyy[ijk]   += Kyy_bh[ijk];
    kyz[ijk]   += Kyz_bh[ijk];
    kzz[ijk]   += Kzz_bh[ijk];

    if(subMinkowsky) {
      alp[ijk] -= 1.;
      gxx[ijk] -= 1.;
      gyy[ijk] -= 1.;
      gzz[ijk] -= 1.;
    }
  }

  // Free memory
  free(psi4_bh);
  free(dpsi4x_bh);
  free(dpsi4y_bh);
  free(dpsi4z_bh);
  free(alpha_bh);
  free(dalphax_bh);
  free(dalphay_bh);
  free(dalphaz_bh);
  free(betax_bh);
  free(betay_bh);
  free(betaz_bh);
  free(gxx_bh);
  free(gxy_bh);
  free(gxz_bh);
  free(gyy_bh);
  free(gyz_bh);
  free(gzz_bh);
  free(Kxx_bh);
  free(Kxy_bh);
  free(Kxz_bh);
  free(Kyy_bh);
  free(Kyz_bh);
  free(Kzz_bh);
}
