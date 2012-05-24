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

//////////////////////////////////////////////////////
//EXTERNAL STANDARD LIBRARIES                         
//////////////////////////////////////////////////////
#ifndef __HDW__
#define __HDW__

//////////////////////////////////////////////////////
//GSL LIBRARIES
//////////////////////////////////////////////////////
#include<gsl/gsl_math.h>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_odeiv.h>

//////////////////////////////////////////////////////
//MACROS
//////////////////////////////////////////////////////
#define PowInt gsl_pow_int
#define Min GSL_MIN
#define Max GSL_MAX

//////////////////////////////////////////////////////
//CUSTOM DATA TYPES
//////////////////////////////////////////////////////

//////////////////////////////////////////////////////
//FUNCTION PROTOTYPE
//////////////////////////////////////////////////////
int HDW_System(double t,const double y[],double dydt[],void *param);

//////////////////////////////////////////////////////
//COMMON CODE
//////////////////////////////////////////////////////
using namespace std;

#endif
