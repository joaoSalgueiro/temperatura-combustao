#include "init.h"

void initAir(CHEMICAL*O2Air,CHEMICAL*N2Air,int P){
	O2Air->got = 0;
	N2Air->got = 0;
	
	O2Air->atoms.C = 0;
	O2Air->atoms.H = 0;
	O2Air->atoms.O = 2;
	O2Air->atoms.N = 0;
	O2Air->atoms.S = 0;
	O2Air->Name = (char*)malloc(3*sizeof(char));
	if(!O2Air->Name){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	strcpy(O2Air->Name,"O2");
	O2Air->Name[2] = '\0';
	O2Air->Mol = 1;
	O2Air->Hf0 = 0;
	
	N2Air->atoms.C = 0;
	N2Air->atoms.H = 0;
	N2Air->atoms.O = 0;
	N2Air->atoms.N = 2;
	N2Air->atoms.S = 0;
	N2Air->Name = (char*)malloc(3*sizeof(char));
	if(!N2Air->Name){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	strcpy(N2Air->Name,"N2");
	N2Air->Name[2] = '\0';
	N2Air->Mol = P*3.76;
	N2Air->Hf0 = 0;
	
	O2Air->cpHighRange.a1 = -1.037939022E+06;
	O2Air->cpHighRange.a2 = 2.344830282E+03;
	O2Air->cpHighRange.a3 = 1.819732036E+00;
	O2Air->cpHighRange.a4 = 1.267847582E-03;
	O2Air->cpHighRange.a5 = -2.188067988E-07;
	O2Air->cpHighRange.a6 = 2.053719572E-11;
	O2Air->cpHighRange.a7 = -8.193467050E-16;
	O2Air->cpHighRange.b1 = -1.689010929E+04;
	
	O2Air->cpLowRange.a1 = -3.425563420E+04;
	O2Air->cpLowRange.a2 =  4.847000970E+02;
	O2Air->cpLowRange.a3 = 1.119010961E+00;
	O2Air->cpLowRange.a4 = 4.293889240E-03;
	O2Air->cpLowRange.a5 = -6.836300520E-07;
	O2Air->cpLowRange.a6 = -2.023372700E-09;
	O2Air->cpLowRange.a7 = 1.039040018E-12;
	O2Air->cpLowRange.b1 = -3.391454870E+03;
	
	N2Air->cpHighRange.a1 = 5.877124060E+05;
	N2Air->cpHighRange.a2 =  -2.239249073E+03;
	N2Air->cpHighRange.a3 = 6.066949220E+00;
	N2Air->cpHighRange.a4 = -6.139685500E-04;
	N2Air->cpHighRange.a5 = 1.491806679E-07;
	N2Air->cpHighRange.a6 = -1.923105485E-11;
	N2Air->cpHighRange.a7 = 1.061954386E-15;
	N2Air->cpHighRange.b1 = 1.283210415E+04;
	
	N2Air->cpLowRange.a1 = 2.210371497E+04;
	N2Air->cpLowRange.a2 = -3.818461820E+02;
	N2Air->cpLowRange.a3 = 6.082738360E+00;
	N2Air->cpLowRange.a4 = -8.530914410E-03;
	N2Air->cpLowRange.a5 = 1.384646189E-05;
	N2Air->cpLowRange.a6 = -9.625793620E-09;
	N2Air->cpLowRange.a7 = 2.519705809E-12;
	N2Air->cpLowRange.b1 = 7.108460860E+02;
}

bool init(CReLIST*p){
	p->count = 0;
	p->maxSize = 12;
	p->List = (CHEMICAL*)malloc(p->maxSize*sizeof(CHEMICAL));
	if(!p->List){
		return True;
	}
	return False;
}

void init_Matrix_of_Functions(){
	/*ptr_Function[0][0] = f0;
	ptr_Function[0][1] = f1;
	ptr_Function[0][2] = f2;
	ptr_Function[0][3] = f3;
	ptr_Function[0][4] = f4;
	ptr_Function[0][5] = f5;
	ptr_Function[0][6] = f6;
	ptr_Function[0][7] = fe1;
	ptr_Function[0][8] = fe2;
	ptr_Function[0][9] = fe3;
	ptr_Function[0][10] = fe4;
	ptr_Function[1][0] = g0;
	ptr_Function[1][1] = g1;
	ptr_Function[1][2] = g2;
	ptr_Function[1][3] = g3;
	ptr_Function[1][4] = g4;
	ptr_Function[1][5] = ge1;
	ptr_Function[1][6] = ge2;
	ptr_Function[1][7] = ge3;
	ptr_Function[2][0] = h0;
	ptr_Function[2][1] = h1;
	ptr_Function[2][2] = h2;
	ptr_Function[2][3] = h3;
	ptr_Function[2][4] = he1;
	ptr_Function[2][5] = he2;
	ptr_Function[2][6] = he3;
	ptr_Function[3][0] = i0;
	ptr_Function[3][1] = i1;
	ptr_Function[3][2] = i2;
	ptr_Function[3][3] = i3;
	ptr_Function[3][4] = i4;
	ptr_Function[3][5] = i5;
	ptr_Function[3][6] = ie1;
	ptr_Function[3][7] = ie2;
	ptr_Function[3][8] = ie3;
	ptr_Function[4][0] = j0a;
	ptr_Function[4][1] = j1a;
	ptr_Function[4][2] = je1;
	ptr_Function[4][3] = je2;*/
	ptr_Function[5][0] = k0;
	ptr_Function[5][1] = k1;
	ptr_Function[5][2] = k2;
	ptr_Function[5][3] = k3;
	ptr_Function[5][4] = ke1;
	ptr_Function[5][5] = ke2;
}
