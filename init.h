#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "newtonMatrix2.h"
#ifndef INIT_H
#define INIT_H

void initAir(CHEMICAL*O2Air,CHEMICAL*N2Air,int P);
bool init(CReLIST*p);
void init_Matrix_of_Functions();

#endif
