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
  parameters[2]=6000;//zc(meters)

  //INITIAL VALUES
  y[0]=Topt;//T
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
  gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(4);
  gsl_odeiv_control *c=gsl_odeiv_control_y_new(0.0,1.0E-7);
  
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
  int it,is;
  double tt;
  double dt=1;//year
  double ddt=dt/100.0;//year
  double Intime=50000.0; //year
  FILE *fl=fopen("HDW-dynstep.dat","w");

  //Basic statistics
  double Tsmin=1000,Tsmed,Tsmax=-1000;
  double awmin=1000,awmed,awmax=-1000;
  int iit;
 
  it=0;
  is=0;
  for(t=0;t<=Intime;t+=dt){

    //if(it>5000){
    if(it>-1){
      fprintf(fl,"%e %e %e %e %e\n",t,y[0],y[1],y[2],y[3]);
    }

    iit=0;
    while(tt<t+dt){
      gsl_odeiv_evolve_apply(e,c,s,&sys,&tt,t+dt,&ddt,y);
      iit++;
    }

    if((it%2000)==0)
      printf("Time t = %e (intermediate steps: %d)\n",t,iit);
    it++;

    if(it>5000){
      Tsmin=Min(y[0],Tsmin);
      Tsmax=Max(y[0],Tsmax);
      Tsmed+=y[0];

      awmin=Min(y[2],awmin);
      awmax=Max(y[2],awmax);
      awmed+=y[2];
      is++;
    }
  }

  Tsmed/=is;
  awmed/=is;
  printf("Tsmin = %e\n",Tsmin-273);
  printf("Tsmed = %e\n",Tsmed-273);
  printf("Tsmax = %e\n",Tsmax-273);

  printf("awmin = %e\n",awmin);
  printf("awmed = %e\n",awmed);
  printf("awmax = %e\n",awmax);

  fclose(fl);
}

/*
  printf("t = %e\n",t);
  printf("\tT = %e\n",y[0]);
  printf("\tac = %e\n",y[1]);
  printf("\taw = %e\n",y[2]);
  printf("\tab = %e\n",y[3]);
*/
