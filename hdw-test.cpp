#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

////////////////////////////////////////////////////
//CONFIGURATION FILE
////////////////////////////////////////////////////
#include <hdw.h>
#include <conf.h>

int HydroDWSystem(double t,const double y[],const double dydt[],void *param);

int main(int argc,char *argv[])
{
  printf("S = %e\n",S);
}
