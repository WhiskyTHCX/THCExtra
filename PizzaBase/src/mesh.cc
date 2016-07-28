#include "mesh.h"

using std::ostream;
using std::endl;


namespace Pizza {

bool region::contains(const gpos &g) const {
  for (int d=0;d<DIMS_GRID;d++) if (g(d)<cmin(d) || g(d)>cmax(d)) return false;
  return true;
}

int region::size() const {
  int s=1;
  for (int k=0;k<DIMS_GRID;k++) s*=size(k);
  return s;
}

mesh::topology mesh::topo;
pz_real mesh::sym_vol;
mesh::extent mesh::ext[3];

void mesh::init(gpos nc,gpos ng0,gpos ng1,vec_u x0,vec_u x1,
  bool obnd0[DIMS_GRID],bool obnd1[DIMS_GRID])
{
  x_min=x0;
  x_max=x1;

  gpos sz = nc+ng0+ng1;

  rall.cmin = 0;
  rall.cmax = sz - 1;
  rint.cmin = ng0;
  rint.cmax = rall.cmax - ng1;
  cindex.init(sz);

  for (int k=0;k<DIMS_GRID;k++) {
    dx(k)             = (x1(k) - x0(k)) / pz_real(nc(k));
    is_outer_bnd0[k]  = obnd0[k];
    is_outer_bnd1[k]  = obnd1[k];
  }
}

pz_real mesh::sym_factor(int dim, int i, bool odd) const
{
  const extent ex = coord_sym(dim);
  if (odd) {
    if (ex==full) return dw(dim);
    return 0;
  }
  if (ex==half) {
    const pz_real p = coord(i,dim)/dw(dim);
    if (p>=0.5) return 2.0*dw(dim);
    if (p<=-0.5) return 0;
    return 2.0*(p+0.5)*dw(dim);
  }
  else if (ex==full) {
    return dw(dim);
  }
  return 1;
}

pz_real mesh::cell_weight(gpos p, bool odd0, bool odd1, bool odd2) const
{
  return sym_vol*sym_factor(0,p(0),odd0)*sym_factor(1,p(1),odd1)*sym_factor(2,p(2),odd2);
}

void mesh::set_topology(topology t, extent e[DIMS_GRID])
{
  topo  = t;
  for (int k=0;k<DIMS_GRID;k++) ext[k] = e[k];
  sym_vol = 1.0;
  if (topo==spherical) {
    error::unless(ext[0]==full, "Coordinate symmetries: combination not implemented.");
    if (ext[2]==flat) sym_vol *= 2;
    if (ext[1]==flat) sym_vol *= (2*M_PI);
  }
  if (topo==cylindrical) {
    error::unless(ext[0]==full && ext[2]!=flat, "Coordinate symmetries: combination not implemented.");
    if (ext[1]==flat) sym_vol *= (2*M_PI);
  }
  if (topo==cartesian) {
    error::unless(ext[0]!=flat && ext[1]!=flat && ext[2]!=flat, "Coordinate symmetries: combination not implemented.");
  }
}

void mesh::coord_spherical(const vec_u& p, pz_real& r, pz_real& th, pz_real& phi) const
{
  switch (coord_type()) {
    case spherical:
      error::unless(coord_sym(0)==full, "PizzaBase: unsupported topology");
      r   = p(0);
      th  = (coord_sym(2)==flat) ? (M_PI / 2.0) : p(2);
      phi = (coord_sym(1)==flat) ? 0 : p(1);
      break;
    case cylindrical:
      error::unless((coord_sym(0)==full) && (coord_sym(2)!=flat),
        "PizzaBase: unsupported topology");
      r    = sqrt(sqr(p(0))+ sqr(p(2)));
      th   = ((p(0) != 0) || (p(2) != 0)) ? atan2(p(0), p(2)) : 0.0;
      phi  = (coord_sym(1)==flat) ? 0 : p(1);
      break;
    case cartesian:
      {
        error::unless((coord_sym(0)!=flat) && (coord_sym(1)!=flat) && (coord_sym(2)!=flat),
          "PizzaBase: unsupported topology");
        const pz_real d    = sqrt(sqr(p(0))+sqr(p(1)));
        r               = p.norm();
        th              = ((d != 0) || (p(2) != 0)) ? atan2(d, p(2)) : 0.0;
        phi             = ((p(0) != 0) || (p(1) != 0)) ? atan2(p(1), p(0)) : 0.0;
      }
      break;
    default:
      assert(false);
  }
}

ostream& operator<<(ostream &o, const mesh &m)
{
  const char* tn[3]={"Cartesian", "Cylindrical", "Spherical"};
  const char* sn[3]={"invariant","mirror","none"};
  o<<"All cells: " << (m.r_all().cmax + 1) << endl
   <<"Interior : " << (m.r_int().cmax - m.r_int().cmin + 1) <<endl
   <<"Ghost0   : "<<m.r_int().cmin <<"  (type [ ";
   for (int d=0;d<DIMS_GRID;d++) o << (m.outerbnd0(d) ? "o " : "i ");
   o<<"])"<<endl
   <<"Ghost1   : "<<(m.r_all().cmax-m.r_int().cmax)<<"  (type [ ";
   for (int d=0;d<DIMS_GRID;d++) o << (m.outerbnd1(d) ? "o " : "i ");
   o<<"])"<<endl
   <<"X_min    : "<<m.xmin()<<endl
   <<"X_max    : "<<m.xmax()<<endl
   <<"Origin   : "<<m.coord(m.r_int().cmin)<<endl
   <<"Spacing  : "<<m.dw()<<endl
   <<"Topology : "<<tn[m.coord_type()] <<endl
   <<"Symmetry : [ ";
   for (int d=0;d<DIMS_GRID;d++) o << sn[m.coord_sym(d)] << " ";
   o<<"]"<<endl;
  return o;
}

ostream& operator<<(ostream &o,const gpos &p)
{
  o<<"["<<p(0);
  for (int k=1;k<DIMS_GRID;k++) o<<","<<p(k);
  return o<<"]";
}

ostream& operator<<(ostream &o,const region &r)
{
  return o<<"{"<<r.cmin<<","<<r.cmax<<"}";
}

}
