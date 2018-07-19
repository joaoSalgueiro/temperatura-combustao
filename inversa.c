#include "inversa.h"
#include "globals.h"

void initializeMatrix(double**matrix,int n){
	int i;
	int j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			matrix[i][j] = 0;
			if(i==j){
				matrix[i][j] = 1;
			}
		}
	}
}

void copy(double**A,double**B,int n){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			B[i][j] = A[i][j];
		}
	}
}

double** inverseMatrix(double**M,int n){
	int i;
	int j;
	int k;
	double pivo;
	double**matrix;
	double**inverted;
	matrix = (double**)malloc(n*sizeof(double*));
	inverted = (double**)malloc(n*sizeof(double*));
	for(i=0;i<n;i++){
		matrix[i] = (double*)malloc(n*sizeof(double));
		inverted[i] = (double*)malloc(n*sizeof(double));
	}
	copy(M,matrix,n);
	initializeMatrix(inverted,n);
	for(j=0;j<n-1;j++){
		pivo = matrix[j][j];
		for(k=0;k<n;k++){
			matrix[j][k]/=pivo;
			inverted[j][k]/=pivo;
		}
		for(i=j+1;i<n;i++){
			pivo = matrix[i][j];
			for(k=0;k<n;k++){
				matrix[i][k]+=matrix[j][k]*(-pivo);
				inverted[i][k]+=inverted[j][k]*(-pivo);
			}
		}
	}
	pivo = matrix[j][j];
	for(k=0;k<n;k++){
		matrix[j][k]/=pivo;
		inverted[j][k]/=pivo;
	}
	for(j=n-1;j>0;j--){
		for(i=0;i<j;i++){
			pivo = matrix[i][j]/matrix[j][j];
			for(k=0;k<n;k++){
				matrix[i][k]+=matrix[j][k]*(-pivo);
				inverted[i][k]+=inverted[j][k]*(-pivo);
			}
		}
	}
	return inverted;
}
