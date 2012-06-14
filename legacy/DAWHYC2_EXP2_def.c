#include <util.c>

int main(int argc,char *argv[])
{
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //VARIABLE DEFINITION
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  float cp,S,sigma,gamma,pp,Topt;
  float q,Ab,Aw,Agf,grad_T,fi,mm,nn,Pmax,EVPmax;
  float Ec,Es,deltat,t_integracion,t_modelacion;
  float Lini,Lfin,deltaL;
  int niter,niterL,nzc;
  Vector zc_vector;
  int zczc,ii;
  float Ac,L,zc;
  float aclouds,awhite,ablack,Ts,d;
  float t,xx,As;
  float k1,k2,k3,k4;
  int j;
  float Tc,Tlb,Tlw,Tanual,I,a,EVP,E;
  float P;
  float k1c,k2c,k3c,k4c;
  float k1w,k2w,k3w,k4w;
  float k1b,k2b,k3b,k4b;
  float Bw,Bb;
  int it;
  Vector time,temperature,white_temperature,
    black_temperature,white_area,black_area,clouds_area,evap,prec;
  Matrix resultados;
  FILE *fl;
  char fname[1000];

  float k1p1,k1p2,k1p3,k1p4,k1p5;

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //PARAMETERS
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cp = 3.e13;	//erg*cm-2*K-1
  S = 2.89e13;	//erg*cm-2*año-1
  sigma = 1789.;	//erg*cm-2*año-1*K-4
  gamma = 0.3;	//año-1
  pp = 1.;	//adimensional
  Topt = 295.5;	//K
  q = 20.;	//K
  Ab = 0.25;	//adimensinal
  Aw = 0.75;	//adimensinal
  Agf = 0.50;	//adimensinal
  grad_T = (-0.0065);	//K*m-1, Trenberth 95, p.10
  
  fi = 0.1;
  mm = 0.35;
  nn = 0.1;
  Pmax = pow((1./mm),(1/nn)); //=36251 [mm/año], cuando ac=1.0
  EVPmax = Pmax;
  //EVPmax = 1511.; //[mm/año] corresponde a Ts=26 grados centígrados
  //Pmax = EVPmax;	//[mm/año]
  //ac_crit = mm*Pmax^nn;

  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //experimento otro: cambiando estos parámetros
  Ec=1.;	//emisividad de las nubes
  Es=1.;	//emisividad de la superficie
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  deltat = 0.01; //[años] tamaño de paso temporal para dTsdt y dadt
  t_integracion = 1; //[año] cada cuanto se guardan valores de las variables
  t_modelacion = 10000;	//[años] período de modelación por cada L
  niter = t_integracion/deltat;	//# de iteraciones en los RK4

  /*
    Lini = 1.992;
    Lfin = 4.004;
    deltaL = 0.004;
  */
  Lini = 1.000;
  Lfin = 1.000;
  deltaL = 0.004;
  niterL = (Lfin-Lini)/deltaL;

  ////zc_vector = [1000.,2000.,3000.,4000.,5000.,6000.,7000.,8000.]
  zc_vector=VecAlloc(1);
  VecSet(zc_vector,0,6000.);
  nzc=1;
    
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  //LOOP IN HEIGHTS
  //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  FILE *ft=fopen("hdw-legacy.dat","w");
  for(zczc=0;zczc<=nzc-1;zczc++){
    zc=VecGet(zc_vector,zczc);
    //Ac=1.-fi*(zc/1000.);
    Ac=0.6;
    resultados=MatAlloc(niterL+1,19);

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    //LOOP IN SOLAR FORCING
    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    for(ii=0;ii<=niterL;ii++){
      L = deltaL*ii+Lini;
      printf("%d de %d, L = %.4lf\n",ii,niterL,L);
      
      //valores iniciales
      aclouds = 0.01	;//adimensional, 0 para reproducir modelo original, 0.01 para iniciar con area de nubes
      awhite = 0.01	;//adimensional
      ablack = 0.01	;//adimensional
      Ts=295.5	;//temperatura en la superficie, valor inicial para rk4

      d = t_modelacion ;//numero de años en el eje de las abscisas - iteraciones de t - dimension de los vectores de resultados

      //printf("Tam:%d\n",(int)(d+1)/2);
      time=VecAlloc((d+1)/2);
      temperature=VecAlloc((d+1)/2);
      white_temperature=VecAlloc((d+1)/2);
      black_temperature=VecAlloc((d+1)/2);
      white_area=VecAlloc((d+1)/2);
      black_area=VecAlloc((d+1)/2);
      clouds_area=VecAlloc((d+1)/2);
      evap=VecAlloc((d+1)/2);
      prec=VecAlloc((d+1)/2);

      it=0;
      
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //LOOP IN MODELLING TIME
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      for(t=0;t<=t_modelacion;t+=t_integracion){

	//if(it>5000){
	if(it>-1){
	  fprintf(ft,"%e %e %e %e %e\n",
		  t,Ts,aclouds,awhite,ablack);

	}

	xx = pp - awhite - ablack;
	As = xx*Agf + ablack*Ab + awhite*Aw;

	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	//TIME INTEGRATION
	//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	for(j=1;j<=niter;j++){
	  
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //TEMPERATURE
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  k1=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*gsl_pow_int((Ts+grad_T*zc),4)-sigma*Es*gsl_pow_int(Ts,4));
	  k2=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*gsl_pow_int(((Ts+k1/2)+grad_T*zc),4)-sigma*Es*gsl_pow_int((Ts+k1/2),4));
	  k3=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*gsl_pow_int(((Ts+k2/2)+grad_T*zc),4)-sigma*Es*gsl_pow_int((Ts+k2/2),4));
	  k4=(1/cp)*(S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As)+sigma*Ec*aclouds*gsl_pow_int(((Ts+k3)+grad_T*zc),4)-sigma*Es*gsl_pow_int((Ts+k3),4));
	  Ts = Ts+deltat*(k1/6+k2/3+k3/3+k4/6);

	  //CLOUD TEMPERATURE
	  Tc=Ts+zc*grad_T;
	  Tlb=q*(As-Ab)+Ts;
	  Tlw=q*(As-Aw)+Ts;
	  
	  //EVAPORATION
	  if(Ts>277){
	    Tanual = Ts - 273.	;//(°C)
	    I = 12.*pow((Tanual/5.),1.5);
	    a = (6.7e-7)*gsl_pow_int(I,3) - (7.7e-5)*PowInt(I,2) + (1.8e-2)*I + 0.49;
	    EVP = 12.*16*pow((10.*(Ts - 273.)/I),a);
	    E = Min(1.,EVP/EVPmax);
	  }else{
	    E = 0.;
	  }

	  //PRECIPITATION
	  P = (1./Pmax)*pow((aclouds/mm),(1./nn));

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //CLOUD COVERING
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  k1c=(1-aclouds)*E-aclouds*P;
	  k2c=(1-(aclouds+k1c/2))*E-(aclouds+k1c/2)*P;
	  k3c=(1-(aclouds+k2c/2))*E-(aclouds+k2c/2)*P;
	  k4c=(1-(aclouds+k3c))*E-(aclouds+k3c)*P;
	  aclouds=aclouds+deltat*(k1c/6+k2c/3+k3c/3+k4c/6);

	  //REPRODUCTIVE FITNESS
	  //WHITE DAISIES
	  if((Tlw>278)&&(Tlw<313))
	    Bw=1-0.003265*PowInt((Topt-Tlw),2);
	  else Bw=0;
	  //BLACK DAISIES
	  if((Tlb>278)&&(Tlb<313)) 
	    Bb=1-0.003265*PowInt((Topt-Tlb),2);
	  else Bb=0;
	  
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //WHITE AREA
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  k1w=awhite*(xx*Bw-gamma);
	  k2w=(awhite+k1w/2)*(xx*Bw-gamma);
	  k3w=(awhite+k2w/2)*(xx*Bw-gamma);
	  k4w=(awhite+k3w)*(xx*Bw-gamma);
	  awhite=awhite+deltat*(k1w/6+k2w/3+k3w/3+k4w/6);

	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  //BLACK AREA
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  k1b=ablack*(xx*Bb-gamma);
	  k2b=(ablack+k1b/2)*(xx*Bb-gamma);
	  k3b=(ablack+k2b/2)*(xx*Bb-gamma);
	  k4b=(ablack+k3b)*(xx*Bb-gamma);
	  ablack=ablack+deltat*(k1b/6+k2b/3+k3b/3+k4b/6);

	  xx = pp - awhite - ablack;
	  As = xx*Agf + ablack*Ab + awhite*Aw;
	  
	}//end for time integration
	
	if(it>5000){
	  VecSet(time,it-5001,t);
	  VecSet(temperature,it-5001,Ts-273);
	  VecSet(white_temperature,it-5001,Tlw-273);
	  VecSet(black_temperature,it-5001,Tlb-273);
	  VecSet(white_area,it-5001,awhite);
	  VecSet(black_area,it-5001,ablack);
	  VecSet(clouds_area,it-5001,aclouds);
	  VecSet(evap,it-5001,E*(1-aclouds));
	  VecSet(prec,it-5001,P*aclouds);
	}
	
	it++;
	
      }//end for modelling time t
      

      if(VERBOSE){
	fprintf(stdout,"Valor %s = %.6e\n",VARS[0],L);
	fprintf(stdout,"Valor %s = %.6e\n",VARS[1],VecMin(temperature))	;//Ts
	fprintf(stdout,"Valor %s = %.6e\n",VARS[2],VecMean(temperature))	;//Ts
	fprintf(stdout,"Valor %s = %.6e\n",VARS[3],VecMax(temperature))	;//Ts
	fprintf(stdout,"Valor %s = %.6e\n",VARS[4],VecMin(white_area))	;//aw
	fprintf(stdout,"Valor %s = %.6e\n",VARS[5],VecMean(white_area))	;//aw
	fprintf(stdout,"Valor %s = %.6e\n",VARS[6],VecMax(white_area))	;//aw
	fprintf(stdout,"Valor %s = %.6e\n",VARS[7],VecMin(black_area))	;//ab
	fprintf(stdout,"Valor %s = %.6e\n",VARS[8],VecMean(black_area))	;//ab
	fprintf(stdout,"Valor %s = %.6e\n",VARS[9],VecMax(black_area))	;//ab
	fprintf(stdout,"Valor %s = %.6e\n",VARS[10],VecMin(clouds_area))	;//ac
	fprintf(stdout,"Valor %s = %.6e\n",VARS[11],VecMean(clouds_area))	;//ac
	fprintf(stdout,"Valor %s = %.6e\n",VARS[12],VecMax(clouds_area))	;//ac
	fprintf(stdout,"Valor %s = %.6e\n",VARS[13],VecMin(evap))	;//E
	fprintf(stdout,"Valor %s = %.6e\n",VARS[14],VecMean(evap))	;//E
	fprintf(stdout,"Valor %s = %.6e\n",VARS[15],VecMax(evap))	;//E
	fprintf(stdout,"Valor %s = %.6e\n",VARS[16],VecMin(prec))	;//P
	fprintf(stdout,"Valor %s = %.6e\n",VARS[17],VecMean(prec))	;//P
	fprintf(stdout,"Valor %s = %.6e\n",VARS[18],VecMax(prec))	;//P
      }
 
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      //RESULTADOS
      //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      MatSet(resultados,ii,0,L);
      MatSet(resultados,ii,1,VecMin(temperature))	;//Ts
      MatSet(resultados,ii,2,VecMean(temperature))	;//Ts
      MatSet(resultados,ii,3,VecMax(temperature))	;//Ts
      MatSet(resultados,ii,4,VecMin(white_area))	;//aw
      MatSet(resultados,ii,5,VecMean(white_area))	;//aw
      MatSet(resultados,ii,6,VecMax(white_area))	;//aw
      MatSet(resultados,ii,7,VecMin(black_area))	;//ab
      MatSet(resultados,ii,8,VecMean(black_area))	;//ab
      MatSet(resultados,ii,9,VecMax(black_area))	;//ab
      MatSet(resultados,ii,10,VecMin(clouds_area))	;//ac
      MatSet(resultados,ii,11,VecMean(clouds_area))	;//ac
      MatSet(resultados,ii,12,VecMax(clouds_area))	;//ac
      MatSet(resultados,ii,13,VecMin(evap))	;//E
      MatSet(resultados,ii,14,VecMean(evap))	;//E
      MatSet(resultados,ii,15,VecMax(evap))	;//E
      MatSet(resultados,ii,16,VecMin(prec))	;//P
      MatSet(resultados,ii,17,VecMean(prec))	;//P
      MatSet(resultados,ii,18,VecMax(prec))	;//P

      VecFree(time);
      VecFree(temperature);
      VecFree(white_temperature);
      VecFree(black_temperature);
      VecFree(white_area);
      VecFree(black_area);
      VecFree(clouds_area);
      VecFree(evap);
      VecFree(prec);

    }//end for ii

    sprintf(fname,"DAWHYC2_EXP2_L0416_%.2f_%d_activa.txt",Ac,(int)zc/1000);
    fl=fopen(fname,"w");
    fprintf(fl,"zc= %.0lf\n",zc);
    fprintf(fl,"Ac= %.2lf\n",Ac);

    for(j=0;j<NVARS;j++)
      fprintf(fl,"%10s ",VARS[j]);

    MatrixFprintf(fl,resultados,"%10.4f ");
    fclose(fl);
    
    MatFree(resultados);

  }//end for heights
  fclose(ft);
  
}//end program
