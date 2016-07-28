#include "eos_barotropic_file.h"
#include "eos_polytropic.h"
#include "eos_tabulated.h"
#include "eos_piecewise_poly.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace Pizza;
using namespace EOS_Barotropic;

template<class T>
T read_param(ifstream& s, const string p)
{
  string line;
  if (!getline(s, line))
    throw runtime_error("EOS Cold: Error while loading from file");

  vector<string> tok;
  boost::split( tok, line, boost::is_any_of("="));
  if ((tok.size() != 2) || (p != boost::trim_copy(tok[0])))
    throw runtime_error("EOS Cold: corrupt datafile format");

  return boost::lexical_cast<T>(boost::trim_copy(tok[1]));
}

/**
Load tabulated EOS from stream. The format is a table of five columns (separated
by whitespace) corresponding to rest mass density \f$ \rho \f$ (SI units),
specific energy \f$\epsilon \f$ (dimensionless),
pressure \f$ P \f$ (SI-units), squared soundspeed  \f$ c_s^2/c^2 \f$,
and \f$ g - 1 \f$ (dimensionless).
See eos_tabulated::eos_tabulated() for more info on the input data.
If there are problems reading the stream or the format is corrupt, an exception is thrown.
*/
eos_1p eos_1p_load_tabulated(ifstream& s, const units& u, const string name)
{
  string l;
  vector<double> v_rmd, v_sed, v_p, v_cs2, v_gm1, v_temp, v_efrac;
  const bool isentropic  = read_param<bool>(s, "isentropic");
  const bool hastemp     = read_param<bool>(s, "has_temp");
  const bool hasefrac    = read_param<bool>(s, "has_efrac");

  while (getline(s,l)) {
    istringstream ls(l);
    double rmd,sed,p,cs2,gm1,temp,efrac;
    ls >> rmd >> sed >> p >> cs2 >> gm1;
    if (hastemp)  ls >> temp;
    if (hasefrac) ls >> efrac;
    if (ls.fail())
      throw runtime_error("EOS Cold: corrupt datafile format");
    v_rmd.push_back(rmd / u.density());
    v_sed.push_back(sed);
    v_p.push_back(p / u.pressure());
    v_cs2.push_back(cs2);
    v_gm1.push_back(gm1);
    if (hastemp)  v_temp.push_back(temp);
    if (hasefrac) v_efrac.push_back(efrac);
  }
  eos_1p eos(new eos_tabulated(v_gm1, v_rmd, v_sed, v_p, v_cs2,
                                  v_temp, v_efrac, isentropic, name));
  return eos;
}

/**
Initialize polytropic EOS from input stream of the format
\verbatim
poly_n    = <n>
poly_rmd  = <rmd>
\endverbatim
where \c \<n\> \c is the polytropic index and \c \<rmd\> \c is the polytropic density scale
given in SI-units. The units to be used by the polytrope can be specified,
but by convention they should satisfy \f$ G=c=1 \f$.
If there are problems reading the stream or the format is corrupt, an exception is thrown.
*/
eos_1p eos_1p_load_polytrope(ifstream& s, const units& u, const string name)
{
  const double poly_n   = read_param<double>(s, "poly_n");
  const double poly_rmd = read_param<double>(s, "poly_rmd") / u.density();
  const double rmd_max  = 1e100; //TODO: do not hardcode polytrope max density
  return eos_1p(new eos_polytrope(poly_n, poly_rmd, rmd_max, name));
}

/**
Initialize piecewise polytropic EOS from input stream of the format
\verbatim
poly_rmd  = <density>
max_rmd   = <density>
<rmd_1>   <gamma_1>
<rmd_2>   <gamma_2>
...
<rmd_n>   <gamma_n>
\endverbatim
where  \c poly_rmd \c is the polytropic density scale of the first
segment and \c max_rmd \c the maximum valid density, both in SI-units.
\c \<rmd_i\> \c are the densities of the lower boundary of each segment.
\c \<gamma_i\> \c are the polytropic exponents of each segment.
The units to be used by the polytrope can be
specified, but by convention they should satisfy \f$ G=c=1 \f$.
If there are problems reading the stream or the format is corrupt,
an exception is thrown.
*/
eos_1p eos_1p_load_piecewise_poly(ifstream& s, const units& u,
                                  const string name)
{
  double poly_rmd = read_param<double>(s, "poly_rmd") / u.density();
  double max_rmd  = read_param<double>(s, "max_rmd") / u.density();
  string l;
  vector<double> v_rmd, v_gamma;
  while (getline(s,l)) {
    istringstream ls(l);
    double rmd, gamma;
    ls >> rmd >> gamma;
    if (ls.fail())
      throw runtime_error("EOS Cold: corrupt datafile format");
    v_rmd.push_back(rmd / u.density());
    v_gamma.push_back(gamma);
  }
  eos_1p eos(new eos_piecewise_poly(poly_rmd,v_rmd, v_gamma,
                                    max_rmd, name));
  return eos;
}


namespace EOS_Barotropic {
/**
Load EOS from a file which starts with the following
\verbatim
name = <name>
type = <type>
\endverbatim
where \c \<name\> \c is a short description of he EOS and
\c \<type\> \c is either polytrope or tabulated, and determines what kind of EOS is returned.
The rest of the file if depends on the type, see eos_1p_load_polytrope()
and eos_1p_load_tabulated().
The units to be used by the resulting EOS can be specified, but should satisfy
\f$ G=c=1 \f$ by convention.
If there are problems loading the file or the format is corrupt, an exception is thrown.
*/
eos_1p load_eos_1p(std::string fname, const units& u)
{
  ifstream f(fname.c_str());
  string name   = read_param<string>(f, "name");
  string type   = read_param<string>(f, "type");
  if (type=="polytrope") {
    return eos_1p_load_polytrope(f, u, name);
  }
  else if (type=="tabulated") {
    return eos_1p_load_tabulated(f, u, name);
  }
  else if (type=="pwpoly") {
    return eos_1p_load_piecewise_poly(f, u, name);
  }
  throw runtime_error("EOS Cold: datafile contains unknown EOS type");
}

}
