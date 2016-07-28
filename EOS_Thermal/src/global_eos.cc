#include "global_eos.h"
#include <stdexcept>

using namespace whizza;

eos_thermal global_eos_thermal::eos;
bool        global_eos_thermal::initialized = false;


eos_thermal global_eos_thermal::get_eos()
{
  if (!initialized) {
    throw std::logic_error("global_eos_thermal: access to eos requested before initialization");
  }
  return eos;
}

void global_eos_thermal::set_eos(const eos_thermal& eos_)
{
  if (initialized) {
    throw std::logic_error("global_eos_thermal: multiple initialization attempted");
  }
  eos         = eos_;
  initialized = true;
}

