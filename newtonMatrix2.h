#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "struct.h"
#include "matrizfunc2.h"

#ifndef NEWTONMATRIX2_H
#define NEWTONMATRIX2_H

double*vectorSub(double*v1,double*v2,int n,double w,int type,int*I);
double**MatrixMult(double**M1,int m1,int n1,double**M2,int m2,int n2);
double*MatrixToVector(double**M,int n);
double*newVector(double*v,double P,double*K,double vC,double vH,double vO,double vN,int size,int type,int*I);
int errorV(double*v1,double*v2,int n,double error,int type);
void getK(double*K,double T,int size,int balance,int*I);
void getInitial(double*Vi,double C,double H,double O,double N,int type);
void completeProducts(CReLIST*p,double*v);
void newProducts(CReLIST*p,double T,double P,double C, double H,double O, double N);

#endif
