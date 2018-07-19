#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "macros.h"
#include "inversa.h"

#ifndef MATRIZFUNC2H
#define MATRIZFUNC2H

double (*ptr_Function[6][11])(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double p,double K,double vC,double vH,double vO,double vN);

double f0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f5(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double f6(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

double g0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double g1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double g2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double g3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double g4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

double h0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double h1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double h2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double h3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

double i0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double i1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double i2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double i3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double i4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double i5(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

double j0a(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double j1a(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

double k0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double k1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double k2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double k3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double ke1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
double ke2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);

void printaMatriz(double**matriz,int n);
double**Jacobian(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double*K,double vC,double vH,double vO,double vN,int size,int type);
double**InverseJacobian(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double*K,double vC,double vH,double vO,double vN,int size,int type);

#endif
