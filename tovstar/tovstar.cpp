/*! \mainpage
The programs tovsol and tovsols are used to compute TOV solutions or solutions
sequences, for arbitrary barotropic equations of state (EOS). For isentropic
EOSs, tovsol can also be used to compute oscillation frequencies and
eigenfunctions in the Cowling approximation (fixed spacetime).
The code is based on a library tovlib. For the latter, python bindings are
available as well for interactive work or scripts.

\par Requirements
  - cmake (make replacement)
  - BOOST libraries
  - GSL library
  - SWIG (optional, needed for Python bindings)
  - Doxygen (optional, needed to build documentation)
  - Python numpy and matplotlib (optional, for plotting of star sequences)
\par Compilation
to compile and install, cd to the source directory and execute
\verbatim
cd build
cmake -DCMAKE_INSTALL_PREFIX:PATH=<install dir> ..
make
make install
\endverbatim
where \<install dir\> is the path where make install installs to.
If the boost library is in a nonstandard location, specify
the -DBOOST_ROOT= argument to cmake.
After installing, make sure \<install dir\>/bin is in your $PATH.
\par Usage of standalone TOV solver
To compute and save a single TOV solution, create a parameter file like
\verbatim
eos         = EOS_BU
rmd_c       = 7.9056e+17
mode_l_max  = 4
mode_n_max  = 6
\endverbatim
where eos is the name of a unified EOS file (
see the tovlib documentation for the format),
rmd_c is the central density in SI units,
and mode_l_max, mode_n_max specifies the maximum order of
the oscillation eigenfunctions to compute.
The order of the lines is fixed and no comments are allowed.
Example EOSs and parameter files can be found in the
\<install dir\>/share/doc/tovstar/examples folder.
To compute the solution, use
\verbatim
tovstar <my_parameter_file>
\endverbatim
This creates a folder star_\<my_parameter_file\>, containing
the following files:
  - description.txt contains general properties of the solution and
    oscillation frequencies.
  - profile.dat contains the star profile
  - the files ef_l*_n*.dat contain the eigenfunctions.

The modes are classified by the spherical harmonic index l and
the number of nodes n. Note n starts at 1 for radial modes, but
0 for nonradial ones, i.e. the lowest radial mode is l=2, n=1,
and the lowest quadrupolar mode is l=2, n=0.

\par Use of standalone mass-radius solver
To compute the mass-radius diagram for a given EOS, create a parameter file of the form
\verbatim
eos         = EOS_BU
min_rmd_c   = 2e16
max_rmd_c   = 8e19
nsteps      = 200
\endverbatim
and then
\verbatim
tovstars <my_parameter_file>
\endverbatim
note the plural s.
The results will be saved in a folder star_seq_\<my_parameter_file\>
as plain text file sequence.dat. You can use the python script
tovstar_plot_seq from within the folder to plot the sequence (use
"tovstar_plot_seq --help" for usage instructions).
*/


#include "modes.h"
#include "eos_barotropic_file.h"
#include <string>
#include <iomanip>
#include <stdexcept>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

namespace bf = boost::filesystem;
namespace po = boost::program_options;
using namespace std;
using namespace Pizza;
using namespace Pizza::TOV;
using namespace EOS_Barotropic;

template<class T>
void getpar(po::variables_map& vm, string name, T& par) {
  if (1 != vm.count(name))
    throw runtime_error(string("You need to specify parameter ")+name);
  par = vm[name].as<T>();
}



class app {
  bool done;
  //double poly_rmd, poly_n;
  double rmd_c, gm1_surf;
  int mode_l_max, mode_n_max;
  bf::path eos_file,odir;
  public:
  app(int argc, char *argv[]);
  void go() const;
};

void app::go() const
{
  if (done) return;
  if (gm1_surf != 0)
    cout << "Warning: using artificial surface cut." <<endl;
  bf::create_directory(odir);

  //eos_cold eos = eos_cold::polytrope(poly_n, poly_rmd);
  eos_1p eos = load_eos_1p(eos_file.string());
  tovsol star(rmd_c, eos, 1.0, 1e-10, gm1_surf);

  bf::path nprof = odir / "profile.dat";
  star.save(nprof.string(), 400, 10.0, units::geom_meter());

  units u = units::si() / units::geom_meter();

  bf:: path ndat = odir / "point.dat";
  bf::ofstream fdat(ndat);
  fdat << setw(18) << (star.rmd_c() / u.density())
       << setw(18) << (star.grav_mass() / u.mass())
       << setw(18) << (star.radius() / u.length())
       << setw(18) << (star.binding_energy() / u.mass())
       << endl;

  bf::path ntxt = odir / "description";
  bf::ofstream ftxt(ntxt);
  if (gm1_surf != 0)
    ftxt << "Warning: artificial stellar surface at gm1 = " << gm1_surf <<endl;
  ftxt << "----- SI-units --------------------" << endl << endl
       << star.to_str(units::si())
       << "----- units G=c=1 -----------------" << endl
       << "(" << units::geom_meter() <<")"<< endl <<endl
       << star.to_str(units::geom_meter()) << endl
       << "----- Frequencies (SI Units) ------" << endl;
  boost::format fmt("l = %d, n = %d, f = %.6e Hz \n");


  for (int l=0; l<=mode_l_max; l++) {
    for (int n= (l>0 ? 0 : 1); n<=mode_n_max; n++) {
      const pmode_cowling mdc(star, l, 0, n, 1e-6);
      const double freq = mdc.freq();
      ftxt << fmt % l % n % (freq / u.freq());

      boost::format md_fn_fmt("ef_l%d_n%d.dat");
      bf::path nef = odir / boost::str(md_fn_fmt % l % n);
      mdc.save_astext(nef.string(), 1000);
    }
  }

}


app::app(int argc, char *argv[]) :
done(false),
odir("star")
{
  po::options_description params("Stellar Parameters");
  params.add_options()
    ("eos", po::value<string>(), "EOS file")
    ("rmd_c", po::value<double>(), "central density / kg m^⁻3")
    ("gm1_surf", po::value<double>()->default_value(0.0), "stop integration when gm1=gm1_surf")
    ("mode_l_max", po::value<int>()->default_value(3), "compute modes up to l=lmax")
    ("mode_n_max", po::value<int>()->default_value(4), "compute modes up to n=nmax");
//    ("poly_rmd", po::value<double>(), "polytropic density scale / kg m^⁻3")
//    ("poly_n", po::value<double>(), "polytropic index")

  po::options_description cmdln("General");
  cmdln.add_options()
    ("help,?", "print help message")
    ("parfile", po::value<string>(), "parameter file");
  cmdln.add(params);

  po::positional_options_description parg;
  parg.add("parfile", 1);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdln).positional(parg).run(), vm);

  if (vm.count("help")) {
    done = true;
    cout
      << "usage: tovstar [options] [parameter file] \n\n"
      << cmdln << endl;
    return;
  }

  if (vm.count("parfile")) {
    bf::path pf(vm["parfile"].as<string>());
    bf::ifstream cfs(pf);
    if (!cfs)
      throw runtime_error("Could not open parameter file");
    po::store(po::parse_config_file(cfs, params), vm);
    odir = bf::path(string("star_")+pf.filename().string());
    eos_file = pf.parent_path();
  }
  notify(vm);

  units u = units::si() / units::geom_meter();
  //getpar(vm, "poly_rmd", poly_rmd);
  //poly_rmd *= u.density();
  //getpar(vm, "poly_n", poly_n);
  string feos;
  getpar(vm, "eos", feos);
  eos_file = eos_file / feos;
  getpar(vm, "rmd_c", rmd_c);
  rmd_c *= u.density();
  getpar(vm, "mode_l_max", mode_l_max);
  getpar(vm, "mode_n_max", mode_n_max);
  getpar(vm, "gm1_surf", gm1_surf);
}


int main(int argc, char *argv[])
{
//  try {
//    try {
      app a(argc, argv);
//      try {
        a.go();
//      }
//      catch (const std::exception& e) {
//        cerr << "Error:" << endl << e.what() << endl;
//      }
//    }
//    catch (const std::exception& e) {
//      cerr << "Parameter error:" << endl << e.what() << endl;
//    }
//  }
//  catch (...) {
//    cerr << "Error: something ugly has happened." << endl;
//  }
}

