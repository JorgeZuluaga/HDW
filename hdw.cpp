/*
#######################################################################
#  _     _ _____   _  _  _ 
# | |   | (____ \ | || || |
# | |__ | |_   \ \| || || |
# |  __)| | |   | | ||_|| |
# | |   | | |__/ /| |___| |
# |_|   |_|_____/  \______|
#
# Hydrological Daisyworld
#######################################################################
# Copyright (C) 2012 Juan F. Salazar, Jorge I. Zuluaga, German Poveda
#######################################################################
*/

#include <hdw.h>
#include <hdw.conf>

int HDW_System(double t,const double y[],double dydt[],void *param)
{
  //FREE PARAMETERES
  double* parameters=(double*)param;
  double L=parameters[0];
  double Ac=parameters[1];
  double zc=parameters[2];
  double As;

  double Ts=y[0];
  double aclouds=y[1];
  double awhite=y[2];
  double ablack=y[3];
  
  double xx;
  xx = 1.0 - awhite - ablack;
  As = xx*Agf + ablack*Ab + awhite*Aw;

  //========================================
  //ENERGY BALANCE
  //========================================
  //y0 = T
  double NSW,ILW,OLW;
  double Tc;

  Tc=Ts+grad_T*zc;
  NSW=S*L*((1-Ac)*aclouds+(1-aclouds))*(1-As);
  ILW=sigma*Ec*aclouds*PowInt(Tc,4);
  OLW=sigma*Es*PowInt(Ts,4);

  dydt[0]=(1/cp)*(NSW+ILW-OLW);

  //========================================
  //CLOUD DYNAMICS
  //========================================
  //y1 = aclouds
  double Tanual,I,a,EVP,E,P;

  //EVAPORATION
  if(Ts>277){
    Tanual = Ts - 273.	;//(°C)
    I = 12.*pow((Tanual/5.),1.5);
    a = (6.7e-7)*PowInt(I,3) - (7.7e-5)*PowInt(I,2) + (1.8e-2)*I + 0.49;
    EVP = 12.*16*pow((10.*(Ts - 273.)/I),a);
    E = Min(1.,EVP/EVPmax);
  }else{
    E = 0.;
  }
  
  //PRECIPITATION
  P = (1./Pmax)*pow((aclouds/mm),(1./nn));

  dydt[1]=(1-aclouds)*E-aclouds*P;

  //========================================
  //DAISY POPULATION DYNAMICS
  //========================================
  //y2 = awhite
  //y3 = ablack
  double Bw,Bb,Tlb,Tlw;

  Tlw=q*(As-Aw)+Ts;
  //GROW RATE
  if((Tlw>278)&&(Tlw<313))
    Bw=1-0.003265*PowInt((Topt-Tlw),2);
  else Bw=0;

  Tlb=q*(As-Ab)+Ts;
  if((Tlb>278)&&(Tlb<313)) 
    Bb=1-0.003265*PowInt((Topt-Tlb),2);
  else Bb=0;

  dydt[2]=awhite*(xx*Bw-Gamma);
  dydt[3]=ablack*(xx*Bb-Gamma);

  return 0;
}
