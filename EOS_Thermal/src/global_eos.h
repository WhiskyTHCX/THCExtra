#ifndef GLOBALEOS_H
#define GLOBALEOS_H

#include "eos_thermal.h"

namespace whizza {

class global_eos_thermal {
  static eos_thermal eos;
  static bool initialized;
  public:
  static eos_thermal get_eos();
  static void set_eos(const eos_thermal& eos_);
};

}

#endif
