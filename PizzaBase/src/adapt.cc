#include "cctk.h"
#include "util_Table.h"
#include "adapt.h"
#include <sstream>
#include <stdexcept>

using std::vector;
using std::ostringstream;
using std::ofstream;
using std::string;

using namespace Pizza;

void cactus_grid::prep(const cGH *cgh) {
  assert(DIMS_GRID==3);
  mcgh=cgh;

  CCTK_INT nbc[6],d1[6],d2[6],d3[6];
  GetBoundarySpecification(6,nbc,d1,d2,d3);

  gpos nc,ghosts0,ghosts1;
  vec_u x0,x1;
  bool outer_bnd0[3];
  bool outer_bnd1[3];

  for (int d=0;d<=2;d++) {
    outer_bnd0[d]  = (cgh->cctk_bbox[2*d  ]==1);
    outer_bnd1[d]  = (cgh->cctk_bbox[2*d+1]==1);
    ghosts0(d)        = outer_bnd0[d] ? nbc[2*d  ] : cgh->cctk_nghostzones[d];
    ghosts1(d)        = outer_bnd1[d] ? nbc[2*d+1] : cgh->cctk_nghostzones[d];
    nc(d)             = cgh->cctk_lsh[d] - ghosts0(d) - ghosts1(d);

    pz_real delta  = cgh->cctk_delta_space[d] / cgh->cctk_levfac[d];
    pz_real origin_glob = cgh->cctk_origin_space[d]
                  + delta * cgh->cctk_levoff[d] / cgh->cctk_levoffdenom[d];

    x0(d) = origin_glob + delta *(pz_real(cgh->cctk_lbnd[d]+ghosts0(d))-0.5);
    x1(d) = x0(d) + delta * pz_real(nc(d));

  }
  refine_factor_ = cgh->cctk_timefac;
  time_          = cgh->cctk_time;
  delta_time_    = cgh->cctk_delta_time / cgh->cctk_timefac;
  iteration_     = cgh->cctk_iteration;

  mesh::init(nc,ghosts0,ghosts1,x0,x1,outer_bnd0,outer_bnd1);
}

pz_real cactus_grid::globalize(pz_real s,const char *method) const {
  pz_real glob;
  int reduction_handle=CCTK_ReductionArrayHandle(method);
  if (reduction_handle<0) CCTK_WARN(0,"Could not get reduction handle");
  int ierr=CCTK_ReduceLocScalar(mcgh,-1,reduction_handle,&s,&glob,CCTK_VARIABLE_REAL);
  if (ierr<0) CCTK_WARN(0,"Could not reduce");
  return glob;
}

pz_cmplx cactus_grid::global_sum(pz_cmplx s) const
{
  const pz_real re = global_sum(s.real());
  const pz_real im = global_sum(s.imag());
  return pz_cmplx(re,im);
}

#include "mpi.h"

int cactus_grid::cpu() const
{
  int c;
  MPI_Comm_rank( MPI_COMM_WORLD, &c );
  return c;
}

void cactus_glop::interpolate(int ivar, const vector<vec_u>& pos,
                               vector<pz_real>& val) const
{
  const size_t npos=pos.size();
  vector<pz_real> coords[3];
  for (size_t d=0; d<3; d++) {
    coords[d].resize(npos);
    for (size_t i=0; i<npos; i++) coords[d][i]=pos[i](d);
  }
  int idop = CCTK_InterpHandle("uniform cartesian");
  error::unless(idop>=0, "Cactus: could not get interpolator handle");
  int idcs = CCTK_CoordSystemHandle("cart3d");
  error::unless(idcs>=0, "Cactus: could not get coord system handle");
  const void* pcoords[3]={&(coords[0][0]),&(coords[1][0]),&(coords[2][0])};
  val.resize(npos);
  void *out[1]={&(val[0])};
  int otyp[1]={CCTK_VARIABLE_REAL};
  int iin[1]={ivar};
  int err = CCTK_InterpGridArrays(mcgh, 3, idop, Util_TableCreateFromString("order=2"),
                    idcs, npos, CCTK_VARIABLE_REAL, pcoords, 1, iin, 1, otyp, out);
  error::unless(err>=0, "Cactus: Interpolate failed");
}

void cactus_glop::interpolate(var_index ivar, const std::vector<vec_u>& pos,
                              std::vector<pz_real>& val) const
{
  interpolate(ivar.index(), pos, val);
}

pz_real cactus_glop::interpolate(int ivar, const vec_u& pos) const
{
  std::vector<vec_u> p(1);
  p[0]=pos;
  std::vector<pz_real> v(1);
  interpolate(ivar,p,v);
  return v[0];
}

pz_real cactus_glop::reduce(const char*op, int ivar) const
{
  pz_real res;
  const int rhandle = CCTK_ReductionHandle(op);
  error::unless(rhandle>=0, "Cactus: could not get reduction handle");
  const int ierr =  CCTK_Reduce (mcgh, -1, rhandle, 1, CCTK_VARIABLE_REAL,
    &res, 1, ivar);
  error::unless(ierr>=0, "Cactus: reduction failed");
  return res;
}

void pzstream::open(string name)
{
  int cpu;
  if (is_open()) return;
  MPI_Comm_rank( MPI_COMM_WORLD, &cpu );
  ostringstream oss;
  oss << name << "_cpu" << cpu;
  clear();
  ofstream::open(oss.str().c_str());
}


var_index::var_index() : idx(-1), ti(0) {}

void var_index::check_valid() const
{
  if (idx<0) {
    throw std::logic_error("Invalid variable index.");
  }
}

void var_index::init(int idx_)
{
  idx   = idx_;
  check_valid();
  ti    = CCTK_VarTypeI(idx);
}

var_index::var_index(std::string name)
{
  int i = CCTK_VarIndex(name.c_str());
  if (i<0) {
    throw std::logic_error(std::string("Cannot find grid function ")
                                        +name);
  }
  init(i);
}

var_index::var_index(int idx_)
{
  init(idx_);
}

var_index::var_index(const var_index& other)
: idx(other.idx), ti(other.ti)
{
  check_valid();
}

const var_index& var_index::operator=(const var_index& other)
{
  idx = other.index();
  ti  = other.type_index();
  return *this;
}

const var_index& var_index::operator=(std::string name)
{
  *this = var_index(name);
  return *this;
}

int var_index::index() const
{
  check_valid();
  return idx;
}

int var_index::type_index() const
{
  check_valid();
  return ti;
}

void var_index::require_type(int ti_) const
{
  if (ti != ti_) {
    throw std::logic_error("Cactus grid function not of required type.");
  }
}

