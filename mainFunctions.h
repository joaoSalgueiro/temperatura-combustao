#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "struct.h"
#include "macros.h"

#ifndef MAINFUNCTIONS_H
#define MAINFUNCTIONS_H

bool Realloc(CReLIST*p);
void ReallocCharMatrix(MATRIX*matrix,int size);
void breakChemical(MATRIX*matrix,char component[]);
ATOM fillComponents(MATRIX matrix);
bool Collect(CReLIST*p,int flag);
double generateTeoricalAir(CReLIST r);
void Combustion(CReLIST r,CReLIST*p,char c,double*air,double*nFuel,double*C,int P);
double formula(Cp cp,double T);
double generateDh(CHEMICAL c,int flag,double n);
int estimateTemperature(CReLIST Reagents,CReLIST Products,CHEMICAL O2Air,CHEMICAL N2Air,double tAir,double nFuel,double*deltaH,int flag);
void getValues(CReLIST*r,CReLIST*p);
void getHf0(CReLIST*p);
double formulaIterativeD(double*v,double T);
double formulaIterative(double*v,double T,double sum);
void printaV(double*v,double sum);
double newtonFormula(double*v,double T, double sum);
void calculateCp(double*v,double temp,CReLIST p);
double newtonIterativeMethod(double temp,double*v,double sum,CReLIST p);
double generateResults(double temp, CReLIST p,double dh);
int saveData(char*name,CReLIST r, CReLIST p,CHEMICAL O2Air,double T,double n,double air);
void unGet(CReLIST*r,CReLIST*p);
void resetProducts(CReLIST*p,double C,double H,double N);

#endif
