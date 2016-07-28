#ifndef ADAPT_H
#define ADAPT_H

#include "mesh.h"
#include "pizza_smatrix.h"
#include <cstdlib>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cstring>
#include <memory>
#include "cGH.h"

namespace Pizza {


template<class T> class pc_lazy {
  T &t;
  size_t i;
  public:
  pc_lazy(T &t_,size_t i_) : t(t_),i(i_) {}
  void operator>>(typename T::base &v) const {t.get(i,v);}
  void operator<<(const typename T::base &v) const {t.set(i,v);}
  void operator<<(const pc_lazy<T> &v) const
    {typename T::base a;v>>a;t.set(i,a);}
};

struct pc_set {
  size_t i;
  pc_set(size_t i_) :i(i_) {}
  template<class T>
  void operator()(T &a,const typename T::base &b) {a.set(i,b);}
};
struct pc_get {
  size_t i;
  pc_get(size_t i_) :i(i_) {}
  template<class T>
  void operator()(T &a,typename T::base &b) {a.get(i,b);}
};
struct pc_fill {
  template<class T>
  void operator()(T &a,const typename T::base &b) {a.fill(b);}
};
struct pc_copy {
  template<class T>
  void operator()(T &a,const T &b) {a.copy(b);}
};

template<class T>
struct  pc_iface : public T {
  void get(size_t i,typename T::base &b) {pc_get g(i); this->apply(g,b);}
  void set(size_t i,const typename T::base &b) {pc_set s(i); this->apply(s,b);}
  void copy(const T &b) {pc_copy c; this->apply(c,b);}
  void fill(const typename T::base &b) {pc_fill f; this->apply(f,b);}
  typedef pc_lazy< pc_iface<T> > etyp;
  etyp operator[](size_t i) {return etyp(*this,i);}
  etyp operator()(const gpos &p) {return etyp(*this,T::m(p));}
};

class scalar_pc {
  pz_real *d;
  map_cart m;
  public:
  typedef pz_real base;
  typedef pz_real& etyp;
  void init(const map_cart &m_,pz_real *d_) {assert(d_); d=d_; m=m_;}
  scalar_pc() :d(0),m() {}
  void get(size_t i,pz_real &b) const {assert(i<m.size()); b=d[i];}
  void set(size_t i,pz_real b) {assert(i<m.size()); d[i]=b;}
  void copy(const scalar_pc &b) const {
    assert(m.size()==b.m.size());
    std::copy(b.d, b.d+m.size(), d);
  }
  void fill(pz_real t) const {std::fill(d,d+m.size(),t);}
  pz_real &operator[](size_t i) {assert(i<m.size()); return d[i];}
  pz_real &operator()(const gpos &p) {return (*this)[m(p)];}
};

class complex_pc {
  pz_real *re;
  pz_real *im;
  map_cart m;
  public:
  typedef pz_cmplx base;
  typedef pc_lazy<complex_pc> etyp;
  void init(const map_cart &m_, pz_real *re_, pz_real* im_) {
    re=re_; im=im_; m=m_;
    assert(re); assert(im);
  }
  complex_pc() :re(0), im(0), m() {}
  void get(size_t i, base& b) const {
    assert(i < m.size()); b = base(re[i], im[i]);
  }
  void set(size_t i, const base& b) {
    assert(i < m.size()); re[i] = b.real(); im[i] = b.imag();
  }
  void copy(const complex_pc &b) const {
    assert(m.size() == b.m.size());
    std::copy(b.re, b.re + m.size(), re);
    std::copy(b.im, b.im + m.size(), im);
  }
  void fill(base t) const {
    std::fill(re, re + m.size(), t.real());
    std::fill(im, im + m.size(), t.imag());
  }
  etyp operator[](size_t i) {return etyp(*this,i);}
  etyp operator()(const gpos &p) {return (*this)[m(p)];}
};

template<int N>
class array_pc {
  pz_real *v[N];
  size_t size;
  public:
  typedef sm_array<N> base;
  array_pc() {size=0; std::fill(v,v+N,(pz_real *)0);}
  void init(size_t size_,pz_real **v_) {
    size=size_;
    for (int k=0;k<N;k++) {v[k]=v_[k]; assert(v[k]);}
  }
  void get(size_t i,base &w) const {
    assert(i<size);
    for (int k=0;k<N;k++) w[k]=v[k][i];
  }
  void set(size_t i,const base &w) const {
    assert(i<size);
    for (int k=0;k<N;k++) v[k][i]=w[k];
  }
  void fill(const base &t) const {
    for (int k=0;k<N;k++)
      std::fill(v[k],v[k]+size,t[k]);
  }
  void copy(const array_pc<N> &b) const {
    assert(size==b.size);
    for (int k=0;k<N;k++)
      std::copy(b.v[k], b.v[k]+size, v[k]);
  }
};

template<int N,bool UP>
struct tensor1_pc_ {
  typedef sm_tensor1<N,UP> base;
  map_cart m;
  array_pc<N> c;
  //tensor1_pc(){}
  void init(const map_cart &m_, pz_real **w) {m=m_; c.init(m_.size(),w);}
  template<class F,class B> void apply(F &f,B &b) {f(c,b.c);}
};

template<bool UP>
struct vec3_pc_ :public tensor1_pc_<3,UP> {
  void init(const map_cart &m, pz_real **w) {tensor1_pc_<3,UP>::init(m,w);}
  void init(const map_cart &m,pz_real *x,pz_real *y,pz_real *z) {
    pz_real *v[]={x,y,z};
    init(m,v);
  }
};

typedef pc_iface<vec3_pc_<true> > vec_u_pc;
typedef pc_iface<vec3_pc_<false> > vec_l_pc;

template<int N,bool UP>
struct tensor2_sym_pc_ {
  map_cart m;
  typedef sm_tensor2_sym<N,UP> base;
  array_pc<base::SIZE> c;
  //tensor2_sym_pc(){}
  void init(const map_cart &m_,pz_real **w) {m=m_; c.init(m_.size(),w);}
  template<class F,class B> void apply(F &f,B &b) {f(c,b.c);}
};

template<bool UP>
struct mats3_pc_ :public tensor2_sym_pc_<3,UP> {
  void init(const map_cart &m,pz_real**w) {tensor2_sym_pc_<3,UP>::init(m,w);}
  void init(const map_cart &m, pz_real*xx, pz_real*xy, pz_real*yy, pz_real*xz, pz_real*yz, pz_real*zz) {
    pz_real*v[]={xx, xy, yy, xz, yz, zz};
    init(m,v);
  }
};

typedef pc_iface<mats3_pc_<true> > mats_u_pc;
typedef pc_iface<mats3_pc_<false> > mats_l_pc;

struct pzstream : public std::ofstream {
  void open(std::string name);
  pzstream() {}
  pzstream(std::string name) {open(name);}
};


class var_index {
  int idx;
  int ti;
  void init(int idx_);
  void check_valid() const;
  public:
  var_index();
  var_index(std::string name);
  explicit var_index(int idx_);
  var_index(const var_index& other);
  const var_index& operator=(const var_index& other);
  const var_index& operator=(std::string name);
  int index() const;
  int type_index() const;
  void require_type(int ti_) const;
};

class cactus_grid : public mesh {
  const cGH *mcgh;
  int iteration_,refine_factor_;
  pz_real time_,delta_time_;
  public:
  void prep(const cGH *cgh);
  pz_real globalize(pz_real s,const char *method) const;
  pz_real global_sum(pz_real s) const {return globalize(s,"sum");}
  pz_real global_min(pz_real s) const {return globalize(s,"minimum");}
  pz_real global_max(pz_real s) const {return globalize(s,"maximum");}
  pz_cmplx global_sum(pz_cmplx s) const;
  int cpu() const;
  bool root() const {return cpu()==0;}
  int refine_factor() const {return refine_factor_;}
  bool coarsest() const {return 1==refine_factor_;}
  pz_real time() const {return time_;}
  pz_real delta_time() const {return delta_time_;}
  int iteration() const {return iteration_;}

};

class cactus_glop {
  const cGH *mcgh;
  public:
  cactus_glop(const cGH *cgh) : mcgh(cgh) {}
  void interpolate(int ivar, const std::vector<vec_u>& pos,
                   std::vector<pz_real>& val) const;
  void interpolate(var_index ivar, const std::vector<vec_u>& pos,
                   std::vector<pz_real>& val) const;
  pz_real interpolate(int ivar, const vec_u& pos) const;
  pz_real reduce(const char*op, int ivar) const;
  pz_real reduce_sum(int ivar) const {return reduce("sum",ivar);}
  pz_real reduce_min(int ivar) const {return reduce("min",ivar);}
  pz_real reduce_max(int ivar) const {return reduce("max",ivar);}
};


template<class T> class cactus_single {
  T* o;
  public:
  cactus_single() : o(0) {}
  ~cactus_single() {if (o) delete o;}
  void operator=(T *o_) {
    if (o) exit(1);
    o=o_;
  }
  T &prep(const cGH *cctkGH) {
    if (0==o) exit(1);
    o->prep(cctkGH);
    return *o;
  }
};

}

#endif
