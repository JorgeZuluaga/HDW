/////////////////////////////////////////////////////
//BEHAVIOR
/////////////////////////////////////////////////////
#define VERBOSE true

/////////////////////////////////////////////////////
//GLOBAL HEADERS
/////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
using namespace std;

/////////////////////////////////////////////////////
//GSL HEADERS
/////////////////////////////////////////////////////
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics_double.h>

/////////////////////////////////////////////////////
//MACROS AND CUSTOM DATATYPES
/////////////////////////////////////////////////////
#define VecAlloc gsl_vector_calloc 
#define MatAlloc gsl_matrix_calloc 

#define VecFree gsl_vector_free
#define MatFree gsl_matrix_free 

#define VecSet gsl_vector_set
#define MatSet gsl_matrix_set
#define VecGet gsl_vector_get
#define MatGet gsl_matrix_get

#define PowInt gsl_pow_int
#define Min GSL_MIN
#define Max GSL_MAX
#define VecMin gsl_vector_min
#define VecMax gsl_vector_max

typedef gsl_vector* Vector;
typedef gsl_matrix* Matrix;

/////////////////////////////////////////////////////
//GLOBAL VARIABLES
/////////////////////////////////////////////////////
#define NVARS 19
char VARS[][1000]={
  "L",
  "Tsmin","Tsmed","Tsmax",
  "awmin","awmed","awmax",
  "abmin","abmed","abmax",
  "acmin","acmed","acmax",
  "Emin","Emed","Emax",
  "Pmin","Pmed","Pmax"
};
/////////////////////////////////////////////////////
//ROUTINES
/////////////////////////////////////////////////////
double VecMean(Vector vec)
{
  return gsl_stats_mean(vec->data,1,(int)vec->size);
}

int MatrixFprintf(FILE *f,Matrix M,const char *frm)
{
  int n=M->size1;
  int m=M->size2;

  //fprintf(stdout,"n=%d,m=%d\n",n,m);
  
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++){
      fprintf(f,frm,MatGet(M,i,j));
    }
    fprintf(f,"\n");
  }
  return 0;
}

/////////////////////////////////////////////////////
//ATTIC
/////////////////////////////////////////////////////
/*
  EXAMPLE CODE VECTOR:

  int *a=NULL;
  a=new int[3];
  
  a[0]=1;

  vector<int> b;
  vector<int> c;
  b.push_back(5);
  c=b;

  printf("SIZE: %d,%d\n",(int)b.size(),c[0]);
  exit(0);
*/
