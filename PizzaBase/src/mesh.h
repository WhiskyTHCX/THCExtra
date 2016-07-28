#ifndef MESH_H
#define MESH_H

#include "pizza_smatrix.h"
#include <cstring>
#include <assert.h>
#include <exception>

namespace Pizza {

enum {DIMS_GRID=3};

struct gpos
{
  int i[DIMS_GRID];

  gpos(int s=0) {for (int k=0;k<DIMS_GRID;k++) i[k]=s;}
  int &operator()(int d) {return i[d];}
  const int &operator()(int d) const {return i[d];}
  gpos &operator+=(const gpos &p)
    {for (int k=0;k<DIMS_GRID;k++) i[k]+=p(k); return *this;}
  gpos &operator-=(const gpos &p)
    {for (int k=0;k<DIMS_GRID;k++) i[k]-=p(k); return *this;}
  bool operator==(const gpos &p1) const;
};

inline gpos operator+(const gpos &a,const gpos &b) {return (gpos(a)+=b);}
inline gpos operator-(const gpos &a,const gpos &b) {return (gpos(a)-=b);}

inline bool gpos::operator==(const gpos &p1) const
{
  for (int k=0;k<DIMS_GRID;k++) if (i[k]!=p1(k)) return false;
  return true;
}

extern std::ostream &operator<<(std::ostream &o,const gpos &p);

struct region {
  struct iterator;

  gpos cmin,cmax;
  region(){}
  region(const gpos &min_,const gpos &max_) : cmin(min_),cmax(max_) {}
  region &range(int dim,int min,int max) {cmin(dim)=min;cmax(dim)=max;return *this;}
  region plane(int dim,int pos) const {return region(*this).range(dim,pos,pos);}
  region face0(int dim) const {return plane(dim,cmin(dim));}
  region face1(int dim) const {return plane(dim,cmax(dim));}
  int size(int dim) const {return cmax(dim)-cmin(dim)+1;}
  int size() const;
  bool operator==(const region &r) const { return (cmin==r.cmin && cmax==r.cmax);}
  bool contains(const gpos &g) const;
};

extern std::ostream &operator<<(std::ostream &o,const region &r);

struct region::iterator: public gpos
{
  const region r;
  iterator(const region &r_) : gpos(r_.cmin),r(r_) {}
  region::iterator &operator++() {
    for (int k=0;(++i[k]>r.cmax(k)) && (k<DIMS_GRID-1);k++) i[k]=r.cmin(k);
    return *this;
  }
  operator bool() const {return i[DIMS_GRID-1]<=r.cmax(DIMS_GRID-1);}
};

class map_cart {
  size_t str[DIMS_GRID];
  public:
  map_cart() {for (int k=0;k<DIMS_GRID;k++) str[k]=0;}
  void init(const gpos &size) {
    str[0]=size(0);
    for (int k=1;k<DIMS_GRID;k++) str[k]=str[k-1]*size(k);
  }
  size_t size() const {return str[DIMS_GRID-1];}
  int operator()(gpos p) const {
    size_t i=p(0);
    for (int k=1;k<DIMS_GRID;k++) i+=str[k-1]*p(k);
    return i;
  }
};

class map_slice {
  int o,s;
  public:
  map_slice(){}
  map_slice(const map_cart &ci, int dim, gpos p) {
    o=ci(p);
    p(dim)++;
    s=ci(p)-o;
  }
  int operator()(int i) {return o + i*s;}
};

class error: public std::exception
{
  const char* msg;
  public:
  error(const char* msg_) : msg(msg_) {}
  virtual const char* what() const throw()
  {
    return msg;
  }
  static void unless(bool c, const char* m)
  {
    if (!c) throw error(m);
  }
  static void incase(bool c, const char* m)
  {
    unless(!c,m);
  }
};

class mesh {
  public:
  enum topology {cartesian, cylindrical, spherical};
  enum extent {flat, half, full};

  private:

  region rall,rint;
  map_cart cindex;
  vec_u x_min,x_max,dx;
  bool is_outer_bnd0[3],is_outer_bnd1[3];
  static topology topo;
  static extent ext[3];
  static pz_real sym_vol;

  pz_real sym_factor(int dim, int i, bool odd) const;

  public:

  const map_cart &mcart() {return cindex;}
  int size() const {return rall.size();}
  void init(gpos nc,gpos ng0,gpos ng1,vec_u x0,vec_u x1,
    bool obnd0[DIMS_GRID],bool obnd1[DIMS_GRID]);

	mesh(){}

  const region &r_all() const {return rall;}
  const region &r_int() const {return rint;}
  region r_ghosts0(int d) const {return region(rint).range(d,0,rint.cmin(d)-1);}
  region r_ghosts1(int d) const {return region(rint).range(d,rint.cmax(d)+1,rall.cmax(d));}

  bool outerbnd0(int dim) const {return is_outer_bnd0[dim];}
  bool outerbnd1(int dim) const {return is_outer_bnd1[dim];}
  int ncells(int dim) const {return rint.size(dim);}
  map_slice slice(const gpos &pos, int dim) const {return map_slice(cindex,dim,pos);}

  const vec_u &dw() const {return dx;}
  pz_real dw(int dim) const {return dx(dim);}
  //pz_real dv() const {pz_real e=1.0; for (int k=0;k<DIMS_GRID;k++) e*=dx(k); return e;}
  const vec_u &xmin() const {return x_min;}
  const vec_u &xmax() const {return x_max;}
  int coord2index(pz_real x, int dim) const
    {return rint.cmin(dim) + int( (x-x_min(dim)) / dx(dim) );}
  gpos coord2index(vec_u x) const {
    gpos e;
    for (int d=0;d<DIMS_GRID;d++) e(d) = coord2index(x(d),d);
    return e;
  }

  pz_real coord(pz_real p,int dim) const {return x_min(dim) + (p-rint.cmin(dim)+0.5)*dx(dim);}
  vec_u coord(gpos p) const
    {vec_u e; for (int k=0;k<DIMS_GRID;k++) e(k)=coord(p(k),k); return e;}
  vec_u origin() const {return coord(rint.cmin);}
  bool samesize(const mesh &m2)
    {return (r_int()==m2.r_int() && r_all()==m2.r_all());}
  static void set_topology(topology t, extent e[DIMS_GRID]);
  static topology coord_type() {return topo;}
  static extent coord_sym(int dim) {return ext[dim];}
  pz_real cell_weight(gpos p, bool odd0=false, bool odd1=false, bool odd2=false) const;
  void coord_spherical(const vec_u& p, pz_real& r, pz_real& th, pz_real& phi) const;
};


extern std::ostream&operator<<(std::ostream &o,const mesh &m);

}

#endif
