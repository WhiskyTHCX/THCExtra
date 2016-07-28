#include <algorithm>
#include <cstring>
#include <fstream>

#include <mpi.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "pizza_central.h"
#include "pizza_idbase.h"
#include "tovsol.h"

#define SQ(X) ((X)*(X))

using namespace Pizza;
using namespace TOV;
using namespace std;

extern "C" void PizzaTOV_BNS_InitialData(CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    bool const set_lapse = CCTK_Equals(initial_lapse, "PizzaTOV");
    bool const set_shift = CCTK_Equals(initial_shift, "PizzaTOV");
    bool const set_temp  = CCTK_Equals(initial_temperature, "PizzaTOV");
    bool const set_efrac = CCTK_Equals(initial_Y_e, "PizzaTOV");

    int const siz = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];

    CCTK_REAL * velx = &vel[0*siz];
    CCTK_REAL * vely = &vel[1*siz];
    CCTK_REAL * velz = &vel[2*siz];

    std::fill(&kxx[0], &kxx[siz], 0.0);
    std::fill(&kxy[0], &kxy[siz], 0.0);
    std::fill(&kxz[0], &kxz[siz], 0.0);
    std::fill(&kyy[0], &kyy[siz], 0.0);
    std::fill(&kyz[0], &kyz[siz], 0.0);
    std::fill(&kzz[0], &kzz[siz], 0.0);
    if(set_shift) {
        std::fill(&betax[0], &betax[siz], 0.0);
        std::fill(&betay[0], &betay[siz], 0.0);
        std::fill(&betaz[0], &betaz[siz], 0.0);
    }

    try {
        units u = Pizza::Base::pizza_base_central::get().internal_units;
        eos_1p eos = Pizza::IDBase::pizza_idbase_central::get().eos;
        tovsol * star[2];
        for(int si = 0; si != 2; ++si) {
            star[si] = new tovsol(star_rhoc[si]/u.density(), eos, 0, 1e-10,
                star_gm1_cut);

            int rank;
            int ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank); assert(0 == ierr);
            if(0 == rank) {
                size_t const ssiz = strlen(out_dir) + strlen("/star[].info") + 2;
                char * fn = new char[ssiz];
                sprintf(fn, "%s/star[%d].info", out_dir, si);
                ofstream os(fn);
                os << star[si]->to_str();
                delete[] fn;
            }
        }

        CCTK_REAL const rs[2] = {
            star[0]->radius(), star[1]->radius()
        };

        CCTK_REAL myx[2], myy[2], myz[2], myr[2];
        bool is_in_star[2];

        CCTK_REAL myrho[2];
        CCTK_REAL myvelx[2], myvely[2], myvelz[2];

        CCTK_REAL mygtt[2];
        CCTK_REAL mygxx[2], mygxy[2], mygxz[2], mygyy[2], mygyz[2], mygzz[2];

        for(int ijk = 0; ijk != siz; ++ijk) {
            /* Solution for each star */
            for(int si = 0; si != 2; ++si) {
                myx[si] = x[ijk] - star_posx[si];
                myy[si] = y[ijk] - star_posy[si];
                myz[si] = z[ijk] - star_posz[si];
                myr[si] = sqrt(SQ(myx[si]) + SQ(myy[si]) + SQ(myz[si]));

                is_in_star[si] = (myr[si] < rs[si]);

                tov_vars const tv = star[si]->operator()(myr[si]);
                CCTK_REAL const fbr2 = 2.0*tv.g_rr()*tv.m_by_r3;

                mygtt[si] = tv.g_tt();

                mygxx[si] = myx[si]*myx[si]*fbr2 + 1;
                mygxy[si] = myx[si]*myy[si]*fbr2;
                mygxz[si] = myx[si]*myz[si]*fbr2;
                mygyy[si] = myy[si]*myy[si]*fbr2 + 1;
                mygyz[si] = myy[si]*myz[si]*fbr2;
                mygzz[si] = myz[si]*myz[si]*fbr2 + 1;

                myrho[si]   = tv.rmd;
                myvelx[si]  = star_velx[si];
                myvely[si]  = star_vely[si];
                myvelz[si]  = star_velz[si];
            }

            /* Combine the two solutions */
            if(set_lapse) {
                CCTK_REAL const gtt = mygtt[0] + mygtt[1] + 1;
                alp[ijk] = sqrt(-gtt);
            }

            gxx[ijk] = mygxx[0] + mygxx[1] - 1;
            gxy[ijk] = mygxy[0] + mygxy[1];
            gxz[ijk] = mygxz[0] + mygxz[1];
            gyy[ijk] = mygyy[0] + mygyy[1] - 1;
            gyz[ijk] = mygyz[0] + mygyz[1];
            gzz[ijk] = mygzz[0] + mygzz[1] - 1;

            CCTK_REAL alpha;
            if(is_in_star[0]) {
                if(is_in_star[1]) {
                    alpha = myrho[0] > myrho[1] ? 1 : 0;
                }
                else {
                    alpha = 1;
                }
            }
            else {
                alpha = 0;
            }
            rho[ijk]  = alpha*myrho [0] + (1 - alpha)*myrho [1];
            velx[ijk] = alpha*myvelx[0] + (1 - alpha)*myvelx[1];
            vely[ijk] = alpha*myvely[0] + (1 - alpha)*myvely[1];
            velz[ijk] = alpha*myvelz[0] + (1 - alpha)*myvelz[1];

            CCTK_REAL const v_x = gxx[ijk]*velx[ijk] + gxy[ijk]*vely[ijk] +
                gxz[ijk]*velz[ijk];
            CCTK_REAL const v_y = gxy[ijk]*velx[ijk] + gyy[ijk]*vely[ijk] +
                gyz[ijk]*velz[ijk];
            CCTK_REAL const v_z = gxz[ijk]*velx[ijk] + gyz[ijk]*vely[ijk] +
                gzz[ijk]*velz[ijk];
            CCTK_REAL const v2 = v_x*velx[ijk] + v_y*vely[ijk] + v_z*velz[ijk];
            w_lorentz[ijk] = sqrt(1.0/(1.0 - v2));

            CCTK_REAL const gm1 = eos.gm1_from_rmd(rho[ijk]);
            eps[ijk] = eos.sed_from_gm1(gm1);
            press[ijk] = eos.p_from_gm1(gm1);
            if(set_temp) {
                temperature[ijk] = eos.temp_from_gm1(gm1);
            }
            if(set_efrac) {
                Y_e[ijk] = eos.efrac_from_gm1(gm1);
            }
        }

        delete star[1];
        delete star[0];
    }
    catch(std::exception & e) {
        CCTK_ERROR(e.what());
    }
}
