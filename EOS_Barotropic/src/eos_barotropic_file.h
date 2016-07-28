#ifndef EOSBAROTROPICFILE_H
#define EOSBAROTROPICFILE_H

#include "eos_barotropic.h"
#include "unitconv.h"
#include <string>


namespace EOS_Barotropic {

///Load EOS from file.
eos_1p load_eos_1p(
  std::string fname,                  ///< Name of file to load from
  const Pizza::units& u=Pizza::units::geom_meter()  ///< Unit system to be used by EOS
);

}

#endif

