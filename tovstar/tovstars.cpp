#include "tovseq.h"
#include "eos_barotropic_file.h"
#include <string>
#include <iomanip>
#include <stdexcept>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#define BOOST_FILESYSTEM_VERSION 3
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

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
  double min_rmd_c, max_rmd_c, gm1_surf;
  int n_steps;
  bf::path eos_file, odir;
  public:
  app(int argc, char *argv[]);
  void go() const;
};

void app::go() const
{
  if (done) return;

  bf::create_directory(odir);

  eos_1p eos = load_eos_1p(eos_file.string());
  tov_sequence sseq(min_rmd_c, max_rmd_c, eos, gm1_surf);
  //tovseq sseq(min_rmd_c, max_rmd_c, eos, n_steps, gm1_surf);

  bf::path nseq       = odir / "sequence.dat";
  bf::path nseq_stbl  = odir / "sequence_stable.dat";
  bf::path nseq_unst  = odir / "sequence_unstable.dat";
  bf::path nstar_maxb = odir / "star_max_bmass";
  bf::path nstar_maxg = odir / "star_max_gmass";

  sseq.save(nseq.string(), nseq_stbl.string(), nseq_unst.string(),
            nstar_maxg.string(), nstar_maxb.string(), n_steps,
            units::si());

}


app::app(int argc, char *argv[]) :
done(false),
odir("star")
{
  po::options_description params("Stellar Parameters");
  params.add_options()
    ("eos", po::value<string>(), "EOS file")
    ("min_rmd_c", po::value<double>(), "minimum central density / kg m^⁻3")
    ("max_rmd_c", po::value<double>(), "maximum central density / kg m^⁻3")
    ("gm1_surf", po::value<double>()->default_value(0.0), "stop integration when hm1=hm1_surf")
    ("nsteps", po::value<int>()->default_value(200), "number of steps to compute");
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
//  getpar(vm, "poly_rmd", poly_rmd);
//  poly_rmd *= u.density();
//  getpar(vm, "poly_n", poly_n);
  string feos;
  getpar(vm, "eos", feos);
  eos_file = eos_file / feos;
  getpar(vm, "min_rmd_c", min_rmd_c);
  min_rmd_c *= u.density();
  getpar(vm, "max_rmd_c", max_rmd_c);
  max_rmd_c *= u.density();
  getpar(vm, "nsteps", n_steps);
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

