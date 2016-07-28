#include "idbase.h"

using namespace Pizza;
using namespace IDBase;


std::auto_ptr<const pizza_idbase_central> pizza_idbase_central::sptr(0);

void pizza_idbase_central::init(const eos_1p& eos_)
{
  error::incase(sptr.get(), "Trying to initialize PizzaBase global variables twice");
  sptr.reset(new pizza_idbase_central(eos_));
}

const pizza_idbase_central& pizza_idbase_central::get()
{
  error::unless(sptr.get(), "Using PizzaBase global variables uninitialized.");
  return *sptr;
}
