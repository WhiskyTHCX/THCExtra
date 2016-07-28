%module pytov
%include "std_string.i"
%include "std_complex.i"
 %{
/* Includes the header in the wrapper code */
#include "eos_barotropic.h"
#include "eos_barotropic_file.h"
#include "eos_polytropic.h"
#include "eos_tabulated.h"
#include "eos_piecewise_poly.h"
#include "tovsol.h"
#include "tovseq.h"
#include "spharmonics.h"
using namespace Pizza;
using namespace Pizza::NumUtils;
using namespace EOS_Barotropic;
 %}
 
 /* Parse the header file to generate wrappers */
 %include "unitconv.h"
 %include "eos_barotropic.h"
 %include "eos_barotropic_file.h"
 %include "eos_polytropic.h"
 %include "eos_tabulated.h"
 %include "eos_piecewise_poly.h"
 %include "tovsol.h"
 %include "tovseq.h"
 %include "spharmonics.h"

 /* ignore eos_cold::load(std::string fname); */

