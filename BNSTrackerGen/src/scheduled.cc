#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "bnstracker.h"
#include "separation.h"

namespace Pizza {
namespace BNSTrackerGen {

void bns_tracker::prep(const cGH *cctkGH)
{
  cactus_grid::prep(cctkGH);
  DECLARE_CCTK_ARGUMENTS;

  const map_cart &m = mcart();
  pz_real *p_glo[]    = {gxx, gxy, gyy, gxz, gyz, gzz};

  rmd_pc.init(m, rho);
  lfac_pc.init(m, w_lorentz);
  glo_pc.init(m, p_glo);
}

cactus_single<bns_tracker> tracker;

}// namespace BNSTrackerGen
}// namespace Pizza



using namespace Pizza;
using namespace Pizza::BNSTrackerGen;
using namespace std;

extern "C" int BNSTrackerGen_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;
  try {
    CCTK_RegisterBanner("BNSTrackerGen: Track binary neutron star positions");
    tracker = new bns_tracker(analysis_reflevel, sym_pi, dist_num_points);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}

extern "C" void BNSTrackerGen_Move_Grids(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (*bns_merger_stage == BNSTrackerGen::UNKNOWN) {
    CCTK_WARN(1, "BNSTrackerGen_MoveGrids: BNS locations not computed yet, skipping.");
    return;
  }
  if (*bns_merger_stage == BNSTrackerGen::INSPIRAL) {
    if (*bns_sep_tot < merge_separation) {
      if (! track_only) {
        active[index_region1]     = 0;
        active[index_region2]     = 0;
        num_levels[index_region3] += add_levels_post_merge;
      }
      *bns_merger_stage         = BNSTrackerGen::MERGING;
    } else {
      if (! track_only) {
        position_x[index_region1] = *bns_x_1;
        position_y[index_region1] = *bns_y_1;
        position_z[index_region1] = 0;

        position_x[index_region2] = *bns_x_2;
        position_y[index_region2] = *bns_y_2;
        position_z[index_region2] = 0;
      }
    }
  }
  if (*bns_merger_stage == BNSTrackerGen::MERGING) {
    if (*bns_sep_tot < collapse_separation) {
      if (! track_only) {
        num_levels[index_region3] += add_levels_post_collapse;
      }

      *bns_merger_stage         = BNSTrackerGen::COLLAPSING;

      *bns_r_1            = 0;
      *bns_x_1            = 0;
      *bns_y_1            = 0;
      *bns_r_2            = 0;
      *bns_x_2            = 0;
      *bns_y_2            = 0;
      *bns_x_md_1         = 0;
      *bns_y_md_1         = 0;
      *bns_r_md_1         = 0;
      *bns_x_md_2         = 0;
      *bns_y_md_2         = 0;
      *bns_r_md_2         = 0;
      *bns_sep_tot        = 0;
      *bns_coord_sep_bc   = 0;
      *bns_proper_sep_bc  = 0;
      *bns_coord_sep_md   = 0;
      *bns_proper_sep_md  = 0;
      *bns_cms_x          = 0;
      *bns_cms_y          = 0;
      *bns_timestamp  = cctk_time;
    }
  }
}

extern "C" void BNSTrackerGen_Track_Stars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if (cctk_levfac[0] != (1 << analysis_reflevel)) return;

  if (((cctk_iteration % analyze_every) != 0) &&
      (*bns_merger_stage != BNSTrackerGen::UNKNOWN)) return;

  if (*bns_merger_stage == BNSTrackerGen::COLLAPSING) return;

  try {
    const bns_locations& locs = tracker.prep(CCTK_PASS_CTOC).track_stars(*bns_phase_tot);

    *bns_phase_1    = locs.star1.phase;
    *bns_r_1        = locs.star1.sep;
    *bns_x_1        = locs.star1.x;
    *bns_y_1        = locs.star1.y;
    *bns_phase_2    = locs.star2.phase;
    *bns_r_2        = locs.star2.sep;
    *bns_x_2        = locs.star2.x;
    *bns_y_2        = locs.star2.y;

    *bns_x_md_1     = locs.star1_maxloc.x;
    *bns_y_md_1     = locs.star1_maxloc.y;
    *bns_r_md_1     = locs.star1_maxloc.sep;
    *bns_phase_md_1 = locs.star1_maxloc.phase;

    *bns_x_md_2     = locs.star2_maxloc.x;
    *bns_y_md_2     = locs.star2_maxloc.y;
    *bns_r_md_2     = locs.star2_maxloc.sep;
    *bns_phase_md_2 = locs.star2_maxloc.phase;

    *bns_phase_tot      = locs.avg.phase;
    *bns_sep_tot        = locs.avg.sep;
    *bns_coord_sep_bc   = locs.coord_sep_bc;
    *bns_coord_sep_md   = locs.coord_sep_md;
    *bns_cms_x          = locs.cms_x;
    *bns_cms_y          = locs.cms_y;
    *bns_timestamp  = cctk_time;

    if (*bns_merger_stage == BNSTrackerGen::UNKNOWN) {
      *bns_merger_stage = BNSTrackerGen::INSPIRAL;
    }
  }
  catch (exception &e) {
    CCTK_WARN(0, e.what());
  }
}

extern "C" void BNSTrackerGen_Init_State(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  *bns_timestamp = -1;
  *bns_merger_stage = BNSTrackerGen::UNKNOWN;
  *bns_phase_tot    = initial_phase;
}

extern "C" void BNSTrackerGen_Separation(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  if ((cctk_iteration % analyze_every) != 0) return;

  if (*bns_merger_stage == BNSTrackerGen::INSPIRAL) {
    cactus_glop glop(CCTK_PASS_CTOC);
    *bns_proper_sep_bc = proper_distance(glop, *bns_x_1, *bns_y_1,
                                   *bns_x_2, *bns_y_2, dist_num_points);
    *bns_proper_sep_md = proper_distance(glop, *bns_x_md_1, *bns_y_md_1,
                                   *bns_x_md_2, *bns_y_md_2, dist_num_points);

  } else {
    *bns_proper_sep_bc = 0;
    *bns_proper_sep_md = 0;
  }
}



