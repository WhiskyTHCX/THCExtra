#ifndef MINIMA_H
#define MINIMA_H

#include "functors.h"

namespace Pizza {
namespace NumUtils {

///Find minumum of a function from initial guess
double find_minimum_from_guess(const function_r2r f, double x0, double x1, double m,
	const double abs_acc, const double rel_acc, int max_iter=1000);

///Find minumum of a function on an interval.
double find_minimum(const function_r2r f, double x0, double x1, const int divisions,
	const double abs_acc, const double rel_acc, int max_iter=1000);

///Find maxumum of a function from initial guess
double find_maximum_from_guess(const function_r2r f, double x0, double x1, double m,
	const double abs_acc, const double rel_acc, int max_iter=1000);

///Find maxumum of a function on an interval.
double find_maximum(const function_r2r f, double x0, double x1, const int divisions,
	const double abs_acc, const double rel_acc, int max_iter=1000);

}
}

#endif
