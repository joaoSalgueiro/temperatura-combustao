#ifndef STRUCTS_H
#define STRUCTS_H

typedef struct _ATOM{
	double C;
	double H;
	double O;
	double N;
	double S;
} ATOM;

typedef struct _Cp{
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double a6;
	double a7;
	double b1;
} Cp;

typedef struct _CHEMICAL{
	ATOM atoms;
	Cp cpLowRange;
	Cp cpHighRange;
	double Mol;
	double Hf0;
	char*Name;
	double Temperature;
	int got;
} CHEMICAL;

typedef struct _CReLIST{
	CHEMICAL*List;
	int count;
	int maxSize;
}CReLIST;

typedef struct _MATRIX{
	char**matrix;
	int count;
} MATRIX;

#endif 
