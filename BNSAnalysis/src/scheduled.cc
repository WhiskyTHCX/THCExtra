#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "bns_analysis.h"

namespace Pizza {
namespace BNSAnalysis {

void bnsanalysis::prep(const cGH *cctkGH)
{
  cactus_grid::prep(cctkGH);
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart();

  size_t s          = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2];
  pz_real *v[]      = {vel, vel+s, vel+2*s};
  pz_real *sc[]      = {scon, scon+s, scon+2*s};
  pz_real *glo[]    = {gxx, gxy, gyy, gxz, gyz, gzz};
  pz_real *beta[]   = {betax, betay, betaz};
  vars_in.init(m, glo, beta, v, sc, w_lorentz, alp, dens, Y_e,
               rho, eps, press);

  if (need_entropy)         p_entropy.init(m, entropy);
  if (set_u_t)              p_u_t.init(m, u_t);
  if (set_s_phi)            p_s_phi.init(m, s_phi);
  if (set_cons_entropy)     p_cons_entropy.init(m, cons_entropy);
  if (set_dens_noatmo)      p_dens_noatmo.init(m, dens_noatmo);
  if (set_dens_unbnd) {
    p_dens_unbnd_geodesic.init(m, dens_unbnd);
    p_dens_unbnd_bernoulli.init(m, dens_unbnd_bernoulli);
    p_dens_unbnd_garching.init(m, dens_unbnd_garching);
  }
  if (set_entropy_unbnd)    p_entropy_unbnd.init(m, entropy_unbnd);
  if (set_ye_unbnd)         p_ye_unbnd.init(m, ye_unbnd);
}

cactus_single<bnsanalysis> bnsanalysis_single;

}
}//namespace Pizza

using namespace Pizza;
using namespace Pizza::BNSAnalysis;
using namespace std;

extern "C" int BNSAnalysis_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;

  try {
    CCTK_RegisterBanner("BNSAnalysis: BNS merger diagnosic quantities");
    double rho_cut = atmo_rho * (1.0 + atmo_tolerance);
    bnsanalysis_single    = new bnsanalysis(set_u_t, set_s_phi,
                                 set_cons_entropy, set_dens_noatmo,
                                 set_dens_unbnd, set_entropy_unbnd,
                                 set_ye_unbnd, rho_cut);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" void BNSAnalysis_Compute(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  if(verbose) {
    CCTK_INFO("BNSAnalysis_Compute");
  }
  try {
    bnsanalysis_single.prep(CCTK_PASS_CTOC).compute();
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

