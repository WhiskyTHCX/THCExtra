#ifndef ROOTS_H
#define ROOTS_H


#include "functors.h"

namespace Pizza {
namespace NumUtils {

///Find root of a function
double findroot(const function_r2r f, double x0, double x1,
	const double abs_acc, const double rel_acc, int max_iter=1000);

}
}

#endif
