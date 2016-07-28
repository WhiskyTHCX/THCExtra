#ifndef CENTRAL_H
#define CENTRAL_H

#include "pizza_unitconv.h"
#include <memory>

namespace Pizza {
namespace Base {

class pizza_base_central
{
  static std::auto_ptr<const pizza_base_central> sptr;
  public:

  units internal_units;
  pizza_base_central(const units& units_)
  : internal_units(units_) {}
  static void init(const units& units_);
  static const pizza_base_central& get();
};


}
}

#endif
