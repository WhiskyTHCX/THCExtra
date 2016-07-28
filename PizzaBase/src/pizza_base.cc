#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "adapt.h"
#include "central.h"
#include <string>

using namespace Pizza;
using std::string;

extern "C" int PizzaBase_Startup(void)
{
  DECLARE_CCTK_PARAMETERS;
  CCTK_RegisterBanner("PizzaBase: pizza infrastructure");

  try {

    mesh::extent ext[3];
    string kwext[3]={coord_sym0, coord_sym1, coord_sym2};
    for (int d=0; d<3; d++) {
      if (kwext[d]=="flat") ext[d]=mesh::flat;
      else if (kwext[d]=="half") ext[d]=mesh::half;
      else ext[d]=mesh::full;
    }

    mesh::topology topo=mesh::cartesian;
    string kwtopo(coord_type);
    if (kwtopo=="spherical") topo=mesh::spherical;
    else if (kwtopo=="cylindrical") topo=mesh::cylindrical;

    mesh::set_topology(topo, ext);

    units glob_units = units::geom_ulength(length_unit);

    Base::pizza_base_central::init(glob_units);
  }
  catch (std::exception &e) {
    CCTK_WARN(0, e.what());
  }
  return 0;
}


