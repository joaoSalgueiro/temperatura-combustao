#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#ifndef INVERSE_H
#define INVERSE_H

void initializeMatrix(double**matrix,int n);
void copy(double**A,double**B,int n);
double** inverseMatrix(double**M,int n);

#endif
