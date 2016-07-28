#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif
int tablereaderc(const char* name);
int tabulated3d(const double indep[3], double dep[], int flag);
double tab3d_press(const double rho, const double eps, const double ye, int* ierr);
double tab3d_csnd2(const double rho, const double eps, const double ye, int* ierr);
double tab3d_temp(const double rho, const double eps, const double ye, int* ierr);
double tab3d_eps(const double rho, const double temp, const double ye);
int eps_range(double myRho, double myYe, double* epsmin, double* epsmax);
void rho_range(double* rho_min, double* rho_max);
void ye_range(double* ye_min, double* ye_max);
void temp_range(double* temp_min, double* temp_max);
double tab3d_entropy_from_eps(const double rho, const double eps, const double ye, int* ierr);
double tab3d_press_from_temp(const double rho, const double temp, const double ye);
double tab3d_csnd2_from_temp(const double rho, const double temp, const double ye);
double tab3d_entropy_from_temp(const double rho, const double temp, const double ye);

#ifdef __cplusplus
}
#endif

int tablereader(const char* name)
{
	return tablereaderc(name);
}


