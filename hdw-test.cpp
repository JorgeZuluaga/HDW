////////////////////////////////////////////////////
//PROGRAM
////////////////////////////////////////////////////
#include <hdw.h>
#include <conf.h>

int main(int argc,char* argv[])
{
  double parameters[3];
  double t,y[4],yerr[4],dydt[4],dydto[4];

  //PARAMETERS
  parameters[0]=1.0;//L(adimensional)
  parameters[1]=0.6;//Ac(adimensional)
  parameters[2]=8000;//zc(meters)

  //INITIAL VALUES
  y[0]=2*Topt;//T
  y[1]=0.01;//ac
  y[2]=0.01;//aw
  y[3]=0.01;//ab

  //DERIVATIVES
  HDW_System(0,y,dydt,parameters);
  printf("dT/dt = %e\n",dydt[0]);
  printf("dac/dt = %e\n",dydt[1]);
  printf("daw/dt = %e\n",dydt[2]);
  printf("dab/dt = %e\n",dydt[3]);

  //LET'S INTEGRATE 
  const gsl_odeiv_step_type *T=gsl_odeiv_step_rk4;
  gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,4);
  gsl_odeiv_system sys;
  sys.function=&HDW_System;
  sys.jacobian=NULL;
  sys.dimension=4;
  sys.params=parameters;

  //INITIAL CODITIONS
  printf("Initial conditions:\n");
  printf("\tT = %e\n",y[0]);
  printf("\tac = %e\n",y[1]);
  printf("\taw = %e\n",y[2]);
  printf("\tab = %e\n",y[3]);

  //t is in years
  double dt=0.1;
  FILE *fl=fopen("HDW.dat","w");
  for(t=0;t<=5000;t+=dt){
    gsl_odeiv_step_apply(s,t,dt,y,yerr,dydt,dydto,&sys);
    dydt[0]=dydto[0];
    dydt[1]=dydto[1];
    dydt[2]=dydto[2];
    dydt[3]=dydto[3];
    /*
    printf("t = %e\n",t);
    printf("\tT = %e\n",y[0]);
    printf("\tac = %e\n",y[1]);
    printf("\taw = %e\n",y[2]);
    printf("\tab = %e\n",y[3]);
    */
    fprintf(fl,"%e %e %e %e %e\n",t,y[0],y[1],y[2],y[3]);
  }
  fclose(fl);
}
