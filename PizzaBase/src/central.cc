#include "central.h"
#include <stdexcept>

using namespace Pizza;
using namespace Pizza::Base;

std::auto_ptr<const pizza_base_central> pizza_base_central::sptr(0);
void pizza_base_central::init(const units& units_)
{
  if (sptr.get()) throw std::runtime_error("Trying to initialize PizzaBase global variables twice");
  sptr.reset(new pizza_base_central(units_));
}

const pizza_base_central& pizza_base_central::get()
{
  if (sptr.get()) return *sptr;
  throw std::runtime_error("Using PizzaBase global variables uninitialized.");
}


