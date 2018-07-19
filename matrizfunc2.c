#include "matrizfunc2.h"

//----------------------------------------------------C + H + N---------------------------------------------------------------------------------

double f0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + 1/(xH2*xH2O*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5));
	}
	return -K + (xH2*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5))/xH2O;
}

double f1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + 1/(pow(xH2,0.5)*(2*vH - xH - 2*xH2 - 2*xH2O)*xH2O*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5));
	}
	return -K + (pow(xH2,0.5)*(2*vH - xH - 2*xH2 - 2*xH2O)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5))/xH2O;
}

double f2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + pow(xO2,0.5)/((vC - vH - 2*vN - xCO2 + xH + 2*xH2 + xH2O + 2*xN2 - xNO2 - 2*xO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5));
	}
	return -K + ((vC - vH - 2*vN - xCO2 + xH + 2*xH2 + xH2O + 2*xN2 - xNO2 - 2*xO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5))/pow(xO2,0.5);
	}

double f3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + pow(xH2,0.5)/(xH*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5));
	}
	return -K + (xH*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5))/pow(xH2,0.5);
}

double f4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + xCO2/((vC - xCO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5));
	}
	return -K + ((vC - xCO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5))/xCO2;
}

double f5(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + (pow(xN2,0.5)*pow(xO2,0.5))/(2*vN - 2*xN2 - xNO2);
	}
	return -K + (2*vN - 2*xN2 - xNO2)/(pow(xN2,0.5)*pow(xO2,0.5));
}

double f6(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	if(K >= 1){
		return -(1/K) + (pow(xN2,0.5)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*xO2)/xNO2;
	}
	return -K + xNO2/(pow(xN2,0.5)*pow(P/(2*vC + vH - xCO2 + xH + xH2 + xN2 - xNO2 - xO2),0.5)*xO2);
}
//----------------------------------------------------C + H + N---------------------------------------------------------------------------------
//------------------------------------------------------C + H-----------------------------------------------------------------------------------

double g0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((xH2*pow(P/(2*vC + vH - xCO2 + xH + xH2 - xO2),0.5)*pow(xO2,0.5))/xH2O,i);
}

double g1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((pow(xH2,0.5)*(2*vH - xH - 2*xH2 - 2*xH2O)*pow(P/(2*vC + vH - xCO2 + xH + xH2 - xO2),0.5))/xH2O,i);
}

double g2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - vH - xCO2 + xH + 2*xH2 + xH2O - 2*xO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 - xO2),0.5))/pow(xO2,0.5),i);
}

double g3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((xH*pow(P/(2*vC + vH - xCO2 + xH + xH2 - xO2),0.5))/pow(xH2,0.5),i);
}

double g4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - xCO2)*pow(P/(2*vC + vH - xCO2 + xH + xH2 - xO2),0.5)*pow(xO2,0.5))/xCO2,i);
}
//------------------------------------------------------C + H-----------------------------------------------------------------------------------
//------------------------------------------------------C + N-----------------------------------------------------------------------------------

double h0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - 2*vN - xCO2 + 2*xN2 - xNO2 - 2*xO2)*pow(P/(2*vC - xCO2 + xN2 - xNO2 - xO2),0.5))/pow(xO2,0.5),i);
	
}

double h1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - xCO2)*pow(P/(2*vC - xCO2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5))/xCO2,i);
}

double h2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((2*vN - 2*xN2 - xNO2)/(pow(xN2,0.5)*pow(xO2,0.5)),i);
}

double h3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(xNO2/(pow(xN2,0.5)*pow(P/(2*vC - xCO2 + xN2 - xNO2 - xO2),0.5)*xO2),i);
}
//------------------------------------------------------C + N-----------------------------------------------------------------------------------
//------------------------------------------------------H + N-----------------------------------------------------------------------------------

double i0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((xH2*pow(P/(vH + xH + xH2 + xN2 - xNO2 - xO2),0.5)*pow(xO2,0.5))/xH2O,i);
}

double i1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((pow(xH2,0.5)*(2*vH - xH - 2*xH2 - 2*xH2O)*pow(P/(vH + xH + xH2 + xN2 - xNO2 - xO2),0.5))/xH2O,i);
}

double i2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((-vH - 2*vN + xH + 2*xH2 + xH2O + 2*xN2 - xNO2 - 2*xO2)*pow(P/(vH + xH + xH2 + xN2 - xNO2 - xO2),0.5))/pow(xO2,0.5),i);
}

double i3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((xH*pow(P/(vH + xH + xH2 + xN2 - xNO2 - xO2),0.5))/pow(xH2,0.5),i);
}

double i4(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow((2*vN - 2*xN2 - xNO2)/(pow(xN2,0.5)*pow(xO2,0.5)),i);
}

double i5(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(xNO2/(pow(xN2,0.5)*pow(P/(vH + xH + xH2 + xN2 - xNO2 - xO2),0.5)*xO2),i);
}
//------------------------------------------------------H + N-----------------------------------------------------------------------------------
//--------------------------------------------------------C-------------------------------------------------------------------------------------

double j0a(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - xCO2 - 2*xO2)*pow(P/(2*vC - xCO2 - xO2),0.5))/pow(xO2,0.5),i);
}

double j1a(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	int i = 1;
	if(K >= 1){
		i = -1;
	}
	return pow(-K,i) + pow(((vC - xCO2)*pow(P/(2*vC - xCO2 - xO2),0.5)*pow(xO2,0.5))/xCO2,i);
}
//--------------------------------------------------------C-------------------------------------------------------------------------------------
//--------------------------------------------------------H-------------------------------------------------------------------------------------

double k0(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("k0: %f\n",-K + (xH2*pow(xO2,0.5)*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/xH2O);
	return -K + (xH2*pow(xO2,0.5)*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/xH2O;
}

double k1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("k1: %f\n",-K + xOH/(pow(xH2,0.5)*pow(xO2,0.5)));
	return -K + xOH/(pow(xH2,0.5)*pow(xO2,0.5));
}

double k2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("k2: %f\n",-K + (xO*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/pow(xO2,0.5));
	return -K + (xO*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/pow(xO2,0.5);
}

double k3(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("k3: %f\n",-K + (xH*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/pow(xH2,0.5));
	return -K + (xH*pow(P/(xH + xH2 + xH2O + xO + xO2 + xOH),0.5))/pow(xH2,0.5);
}

double ke1(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("ke1: %f\n",-vH + xH + 2*xH2 + 2*xH2O + xOH);
	return -vH + xH + 2*xH2 + 2*xH2O + xOH;
}

double ke2(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN){
	//printf("ke2: %f\n",-vO + xH2O + xO + 2*xO2 + xOH);
	return -vO + xH2O + xO + 2*xO2 + xOH;
}
//--------------------------------------------------------H-------------------------------------------------------------------------------------

void printaMatriz(double**matriz,int n){
	int i;
	int j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			printf("|%f|",matriz[i][j]);
		}
		printf("\n");
	}exit(1);
}

double**Jacobian(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double*K,double vC,double vH,double vO,double vN,int size,int type){
	int i,j;
	//puts("1");
	double**J;
	//puts("2");
	
	J = (double**)malloc(size*sizeof(double*));
	for(i = 0;i < size;i++){
		J[i] = (double*)malloc(size*sizeof(double));
	}
	for(i = 0;i < size;i++){
		switch(type){
			/*case 1:
			{
				double (*ptr_Function[size]) (double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);				
				ptr_Function[0] = _f0;
				ptr_Function[1] = _f1;
				ptr_Function[2] = _f2;
				ptr_Function[3] = _f3;
				ptr_Function[4] = f4;
				ptr_Function[5] = f5;
				ptr_Function[6] = f6;
				J[i][0] = ((*ptr_Function[i]) (xCO2 + dx,xH2O,xH2,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2 - dx,xH2O,xH2,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][1] = ((*ptr_Function[i]) (xCO2,xH2O + dx,xH2,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O - dx,xH2,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][2] = ((*ptr_Function[i]) (xCO2,xH2O,xH2 + dx,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O,xH2 - dx,xH,xO2,xN2,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][3] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH + dx,xO2,xN2,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH - dx,xO2,xN2,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][4] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2 + dx,xN2,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2 - dx,xN2,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][5] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2,xN2 + dx,xNO2,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2,xN2 - dx,xNO2,P,K[i],vC,vH,vN))/(2*dx);
				J[i][6] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2,xN2,xNO2 + dx,P,K[i],vC,vH,vN) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2,xN2,xNO2 - dx,P,K[i],vC,vH,vN))/(2*dx);
				break;
			}
			case 2:
			{
				//printf("%d\n",i);
				double (*ptr_Function[size]) (double xCO2,double xH2O,double xH2,double xH,double xO2,double P,double K,double vC,double vH);
				ptr_Function[0] = g0;
				ptr_Function[1] = g1;
				ptr_Function[2] = g2;
				ptr_Function[3] = g3;
				ptr_Function[4] = g4;
				J[i][0] = ((*ptr_Function[i]) (xCO2 + dx,xH2O,xH2,xH,xO2,P,K[i],vC,vH) - (*ptr_Function[i]) (xCO2 - dx,xH2O,xH2,xH,xO2,P,K[i],vC,vH))/(2*dx);
				J[i][1] = ((*ptr_Function[i]) (xCO2,xH2O + dx,xH2,xH,xO2,P,K[i],vC,vH) - (*ptr_Function[i]) (xCO2,xH2O - dx,xH2,xH,xO2,P,K[i],vC,vH))/(2*dx);
				J[i][2] = ((*ptr_Function[i]) (xCO2,xH2O,xH2 + dx,xH,xO2,P,K[i],vC,vH) - (*ptr_Function[i]) (xCO2,xH2O,xH2 - dx,xH,xO2,P,K[i],vC,vH))/(2*dx);
				J[i][3] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH + dx,xO2,P,K[i],vC,vH) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH - dx,xO2,P,K[i],vC,vH))/(2*dx);
				J[i][4] = ((*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2 + dx,P,K[i],vC,vH) - (*ptr_Function[i]) (xCO2,xH2O,xH2,xH,xO2 - dx,P,K[i],vC,vH))/(2*dx);
				break;
			}
			case 3:
			{
				double (*ptr_Function[size]) (double xCO2,double xO2,double xN2,double xNO2,double P,double K,double vC,double vN);
				ptr_Function[0] = h0;
				ptr_Function[1] = h1;
				ptr_Function[2] = h2;
				ptr_Function[3] = h3;
				J[i][0] = ((*ptr_Function[i]) (xCO2 + dx,xO2,xN2,xNO2,P,K[i],vC,vN) - (*ptr_Function[i]) (xCO2 - dx,xO2,xN2,xNO2,P,K[i],vC,vN))/(2*dx);
				J[i][1] = ((*ptr_Function[i]) (xCO2,xO2 + dx,xN2,xNO2,P,K[i],vC,vN) - (*ptr_Function[i]) (xCO2,xO2 - dx,xN2,xNO2,P,K[i],vC,vN))/(2*dx);
				J[i][2] = ((*ptr_Function[i]) (xCO2,xO2,xN2 + dx,xNO2,P,K[i],vC,vN) - (*ptr_Function[i]) (xCO2,xO2,xN2 - dx,xNO2,P,K[i],vC,vN))/(2*dx);
				J[i][3] = ((*ptr_Function[i]) (xCO2,xO2,xN2,xNO2 + dx,P,K[i],vC,vN) - (*ptr_Function[i]) (xCO2,xO2,xN2,xNO2 - dx,P,K[i],vC,vN))/(2*dx);
				break;
			}
			case 4:
			{
				double (*ptr_Function[size]) (double xH2O,double xH2,double xH,double xO2,double xN2,double xNO2,double P,double K,double vH,double vN);
				ptr_Function[0] = i0;
				ptr_Function[1] = i1;
				ptr_Function[2] = i2;
				ptr_Function[3] = i3;
				ptr_Function[4] = i4;
				ptr_Function[5] = i5;
				J[i][0] = ((*ptr_Function[i]) (xH2O + dx,xH2,xH,xO2,xN2,xNO2,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O - dx,xH2,xH,xO2,xN2,xNO2,P,K[i],vH,vN))/(2*dx);
				J[i][1] = ((*ptr_Function[i]) (xH2O,xH2 + dx,xH,xO2,xN2,xNO2,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O,xH2 - dx,xH,xO2,xN2,xNO2,P,K[i],vH,vN))/(2*dx);
				J[i][2] = ((*ptr_Function[i]) (xH2O,xH2,xH + dx,xO2,xN2,xNO2,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O,xH2,xH - dx,xO2,xN2,xNO2,P,K[i],vH,vN))/(2*dx);
				J[i][3] = ((*ptr_Function[i]) (xH2O,xH2,xH,xO2 + dx,xN2,xNO2,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O,xH2,xH,xO2 - dx,xN2,xNO2,P,K[i],vH,vN))/(2*dx);
				J[i][4] = ((*ptr_Function[i]) (xH2O,xH2,xH,xO2,xN2 + dx,xNO2,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O,xH2,xH,xO2,xN2 - dx,xNO2,P,K[i],vH,vN))/(2*dx);
				J[i][5] = ((*ptr_Function[i]) (xH2O,xH2,xH,xO2,xN2,xNO2 + dx,P,K[i],vH,vN) - (*ptr_Function[i]) (xH2O,xH2,xH,xO2,xN2,xNO2 - dx,P,K[i],vH,vN))/(2*dx);
				break;
			}
			case 5:{
				double (*ptr_Function[size]) (double xCO2,double xO2,double P,double K,double vC);
				ptr_Function[0] = J0;
				ptr_Function[1] = J1;
				J[i][0] = ((*ptr_Function[i]) (xCO2 + dx,xO2,P,K[i],vC) - (*ptr_Function[i]) (xCO2 - dx,xO2,P,K[i],vC))/(2*dx);
				J[i][1] = ((*ptr_Function[i]) (xCO2,xO2 + dx,P,K[i],vC) - (*ptr_Function[i]) (xCO2,xO2 - dx,P,K[i],vC))/(2*dx);
				break;
			}*/
			case 6:
			{
				J[i][0] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O + dx,xH2,xH,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O - dx,xH2,xH,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				J[i][1] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2 + dx,xH,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2 - dx,xH,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				J[i][2] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH + dx,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH - dx,xO2,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				J[i][3] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2 + dx,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2 - dx,xO,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				J[i][4] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2,xO + dx,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2,xO - dx,xOH,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				J[i][5] = ((*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2,xO,xOH + dx,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN) - (*ptr_Function[type - 1][i]) (xCO2,xCO,xH2O,xH2,xH,xO2,xO,xOH - dx,xN2,xNO,xNO2,P,K[i],vC,vH,vO,vN))/(2*dx);
				break;
			}
		}					
		
	}
	return J;
}

double**InverseJacobian(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double*K,double vC,double vH,double vO,double vN,int size,int type){
	double**M = Jacobian(xCO2,xCO,xH2O,xH2,xH,xO2,xO,xOH,xN2,xNO,xNO2,P,K,vC,vH,vO,vN,size,type);
	
	//printaMatriz(M,size);
	M = inverseMatrix(M,size);
	return M;
}
