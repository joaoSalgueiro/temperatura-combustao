#include "newtonMatrix2.h"
double*vectorSub(double*v1,double*v2,int n,double w,int type,int*I){
	int i;
	
	double*v = (double*)malloc(max*sizeof(double));
	for(i = 0; i < max;i++){
		v[i] = 0;
	}
	for(i = 0;i < n;i++){
		v[I[i]] = v1[I[i]] - w*v2[i];
	}
	return v;
}

double**MatrixMult(double**M1,int m1,int n1,double**M2,int m2,int n2){
	int i;
	int j;
	int k;
	if(n1 != m2){
		return 0;
	}
	double**M = (double**)malloc(m1*sizeof(double*));
	for(i = 0; i < m1;i++){
		M[i] = (double*)malloc(n2*sizeof(double));
	}
	for(i=0;i<m1;i++){
		for(j=0;j<n2;j++){
			M[i][j] = 0;
		}
	}
	for(k = 0;k<m1;k++){
		for(i = 0;i<n2;i++){
			for(j = 0;j<m2;j++){
				M[k][i] += M1[k][j]*M2[j][i];
			}
		}
	}
	return M;
}

double*MatrixToVector(double**M,int n){
	int i;
	double*V = (double*)malloc(n*sizeof(double));
	for(i = 0;i < n;i++){
		V[i] = M[i][0];
	}
	return V;
}
/*
double*newVector(double*v,double P,double*K,double vC,double vH,double vO,double vN,int size,int type,int*I){
	int i;
	
	double**F = (double**)malloc(size*sizeof(double*));
	for(i = 0;i < size;i++){
		F[i] = (double*)malloc(sizeof(double));
	}
	
	switch(type){
		/*case 1:
		{
			double (*ptr_Function[size]) (double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
			ptr_Function[type][0] = f0;
			ptr_Function[type][1] = f1;
			ptr_Function[type][2] = f2;
			ptr_Function[type][3] = f3;
			ptr_Function[type][4] = f4;
			ptr_Function[type][5] = f5;
			ptr_Function[type][6] = f6;
			for(i = 0;i < size;i++){
				F[i][0] = (*ptr_Function[i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vCO2,vH2O,vN2);
			}
			break;
		}
		case 2:
		{
			double (*ptr_Function[size]) (double xCO2,double xH2O,double xH2,double xH,double xO2,double P,double K,double vCO2,double vH2O);
			ptr_Function[0] = g0;
			ptr_Function[1] = g1;
			ptr_Function[2] = g2;
			ptr_Function[3] = g3;
			ptr_Function[4] = g4;
			for(i = 0;i < size;i++){
				F[i][0] = (*ptr_Function[i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vCO2,vH2O,vN2);
			}
			break;
		}
		case 3:
		{
			double (*ptr_Function[size]) (double xCO2,double xO2,double xN2,double xNO2,double P,double K,double vCO2,double vN2);
			ptr_Function[0] = h0;
			ptr_Function[1] = h1;
			ptr_Function[2] = h2;
			ptr_Function[3] = h3;
			for(i = 0;i < size;i++){
				F[i][0] = (*ptr_Function[i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vCO2,vH2O,vN2);
			}
			break;
		}
		case 4:
		{
			double (*ptr_Function[size]) (double xH2O,double xH2,double xH,double xO2,double xN2,double xNO2,double P,double K,double vH2O,double vN2);
			ptr_Function[0] = i0;
			ptr_Function[1] = i1;
			ptr_Function[2] = i2;
			ptr_Function[3] = i3;
			ptr_Function[4] = i4;
			ptr_Function[5] = i5;
			for(i = 0;i < size;i++){
				F[i][0] = (*ptr_Function[i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vH2O,vN2);
			}
			break;
		}
		case 5:
		{
			double (*ptr_Function[size]) (double xCO2,double xO2,double P,double K,double vCO2);
			ptr_Function[0] = J0;
			ptr_Function[1] = J1;
			

			break;
		}
		case 6:
		{
			double (*ptr_Function[size])(double xCO2,double xCO,double xH2O,double xH2,double xH,double xO2,double xO,double xOH,double xN2,double xNO,double xNO2,double P,double K,double vC,double vH,double vO,double vN);
			ptr_Function[0] = k0;
			ptr_Function[1] = k1;
			ptr_Function[2] = k2;
			ptr_Function[3] = k3;
			ptr_Function[4] = ke1;
			ptr_Function[5] = ke2;
			for(i = 0;i < size;i++){
				F[i][0] = (*ptr_Function[i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vC,vH,vO,vN);
				//printf("F: %f\n",F[i][0]);
			}//getchar();
			break;
		}
	}
	double**J = InverseJacobian(v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K,vC,vH,vO,vN,size,type);
	//puts("aqui^");
	double**M = MatrixMult(J,size,size,F,size,1);
	double*V = MatrixToVector(M,size);
	return vectorSub(v,V,size,relaxing,type,I);
}*/

double*newVector(double*v,double P,double*K,double vC,double vH,double vO,double vN,int size,int type,int*I){
	int i;
	
	double**F = (double**)malloc(size*sizeof(double*));
	for(i = 0;i < size;i++){
		F[i] = (double*)malloc(sizeof(double));
	}
	for(i = 0;i < size;i++){
		F[i][0] = (*ptr_Function[type - 1][i]) (v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K[i],vC,vH,vO,vN);
	}
	double**J = InverseJacobian(v[0],v[7],v[1],v[2],v[3],v[4],v[9],v[10],v[5],v[8],v[6],P,K,vC,vH,vO,vN,size,type);
	//puts("aqui^");
	double**M = MatrixMult(J,size,size,F,size,1);
	double*V = MatrixToVector(M,size);
	return vectorSub(v,V,size,relaxing,type,I);
}

int errorV(double*v1,double*v2,int n,double error,int type){
	int i;
	int*I = (int*)malloc(n*sizeof(int));
	switch(type){
		case 1:
			I[0] = 0;I[1] = 1;I[2] = 2;I[3] = 3;I[4] = 4;I[5] = 5;I[6] = 6;
			break;
		case 2:
			I[0] = 0;I[1] = 1;I[2] = 2;I[3] = 3;I[4] = 4;
			break;
		case 3:
			I[0] = 0;I[1] = 4;I[2] = 5;I[3] = 6;
			break;
		case 4:
			I[0] = 1;I[1] = 2;I[2] = 3;I[3] = 4;I[4] = 5;I[5] = 6;
			break;
		case 5:
			I[0] = 0;I[1] = 4;
			break;
		case 6:
			I[0] = 1;I[1] = 2;I[2] = 3;I[3] = 4;
			break;
	}
	for(i = 0;i < n;i++){
		if(abs(v1[I[i]] - v2[I[i]]) >= error){
			return 1;
		}
	}
	return 0;
}
	
void getK(double*K,double T,int size,int balance,int*I){
	int i;
	FILE*arq = fopen("table4.txt","r");
	double T1, T2, dummy;
	char*line1 = (char*)malloc(100*sizeof(char));;
	char*line2 = (char*)malloc(100*sizeof(char));;
	double*k1 = (double*)malloc(max*sizeof(double));
	double*k2 = (double*)malloc(max*sizeof(double));
	fgets(line1,100,arq);
	fgets(line2,100,arq);
	sscanf(line1,"%lf",&T1);
	sscanf(line2,"%lf",&T2);
	while(T > T2){
		for(i = 0;i < 100;i++){
			line1[i] = line2[i];
		}
		fgets(line2,100,arq);
		sscanf(line1,"%lf",&T1);
		sscanf(line2,"%lf",&T2);
		
	}
	//printf("%f %f\n",T1,T2);exit(1);
	sscanf(line1,"%lf %lf %lf %lf %lf %lf %lf %lf",&dummy,&k1[0],&k1[1],&k1[2],&k1[3],&k1[4],&k1[5],&k1[6]);
	sscanf(line2,"%lf %lf %lf %lf %lf %lf %lf %lf",&dummy,&k2[0],&k2[1],&k2[2],&k2[3],&k2[4],&k2[5],&k2[6]);
	for(i = 0;i < size - balance;i++){
		K[i] = (k2[I[i]] - k1[I[i]])*T/(T2 - T1) + k1[I[i]] - (k2[I[i]] - k1[I[i]])*T1/(T2 - T1);
	}
	for(i = 0;i < size - balance;i++){
		K[i] = pow(10,K[i]);
	}
}
void getInitial(double*Vi,double C,double H,double O,double N,int type){
	switch(type){
		/*case 1:
			Vi[0] = C - iV;Vi[1] = H - iV;Vi[2] = iV;Vi[3] = iV;Vi[4] = iV;Vi[5] = N - iV;Vi[6] = iV;
			break;
		case 2:
			Vi[0] = C - iV;Vi[1] = H - iV;Vi[2] = iV;Vi[3] = iV;Vi[4] = iV;
			break;
		case 3:
			Vi[0] = C - iV;Vi[4] = iV;Vi[5] = N - iV;Vi[6] = iV;
			break;
		case 4:
			Vi[1] = H - iV;Vi[2] = iV;Vi[3] = iV;Vi[4] = iV;Vi[5] = N - iV;Vi[6] = iV;
			break;
		case 5:
			Vi[0] = C - iV;Vi[4] = iV;
			break;*/
		case 6:
			//H = 4;
			Vi[1] = H - iV;Vi[2] = iV;Vi[3] =iV;Vi[4] = O - iV; Vi[9] = iV; Vi[10] = iV;
			break;
	}
}

void completeProducts(CReLIST*p,double*v){
	int i;
	p->List[0].Mol = v[0];
	p->List[1].Mol = v[1];
	p->List[2].Mol = v[2];
	p->List[3].Mol = v[3];
	p->List[4].Mol = v[4];
	p->List[5].Mol = v[5];
	p->List[6].Mol = v[6];
	p->List[8].Mol = v[7];
	p->List[9].Mol = v[8];
	p->List[10].Mol = v[9];
	p->List[11].Mol = v[10];
	for(i = 0;i < 12;i++){
		printf("%f\n",p->List[i].Mol);
	}printf("\n");getchar();
}

void newProducts(CReLIST*p,double T,double P,double C, double H,double O, double N){
	int size;
	int type;
	int balance;
	int sum = 0;
	if (C){ sum += 'C';}
	if (H){ sum += 'H';}
	if (N){ sum += 'N';}
	switch(sum){
		case 217:
			type = 1;
			size = 7;
			break;
		case 139:
			type = 2;
			size = 5;
			break;
		case 145:
			type = 3;
			size = 4;
			break;
		case 150:
			type = 4;
			size = 6;
			break;
		case 67:
			type = 5;
			size = 2;
			break;
		case 72:
			type = 6;
			size = 6;
			balance = 2;
			break;
	}
	int i;
	/*for(i = 0;i < 7;i++){
		printf("%f\n",K[i]);
	}exit(1);*/
	int*I = (int*)malloc(size*sizeof(int));
	int*I2 = (int*)malloc(size*sizeof(int));
	switch(type){
		/*case 1:
			I[0] = 0;I[1] = 1;I[2] = 2;I[3] = 3;I[4] = 4;I[5] = 5;I[6] = 6;
			I2[0] = 2;I2[1] = 3;I2[2] = 1;I2[3] = 0;I2[4] = 4;I2[5] = 5;I2[6] = 6;
			break;
		case 2:
			I[0] = 0;I[1] = 1;I[2] = 2;I[3] = 3;I[4] = 4;
			I2[0] = 2;I2[1] = 3;I2[2] = 1;I2[3] = 0;I2[4] = 4;
			break;
		case 3:
			I[0] = 0;I[1] = 4;I[2] = 5;I[3] = 6;
			I2[0] = 1;I2[1] = 4;I2[2] = 5;I2[3] = 6;
			break;
		case 4:
			I[0] = 1;I[1] = 2;I[2] = 3;I[3] = 4;I[4] = 5;I[5] = 6;
			I2[0] = 2;I2[1] = 3;I2[2] = 1;I2[3] = 0;I2[4] = 5;I2[5] = 6;
			break;
		case 5:
			I[0] = 0;I[1] = 4;
			I2[0] = 1;I2[1] = 4;
			break;*/
		case 6:
			I[0] = 1;I[1] = 2;I[2] = 3;I[3] = 4; I[4] = 9;I[5] = 10;
			I2[0] = 2;I2[1] = 3;I2[2] = 0;I2[3] = 1;
			break;
	}
	if(testing){
		T = Ttest;
		P = Ptest;
	}
	double*K = (double*)malloc(max*sizeof(double));
	for(i = 0;i < max;i++){
		K[i] = 0;
	}
	//puts("ok");
	getK(K,T,size,balance,I2);
	for(i = 0;i < max;i++){
		printf("Kp: %f\n",K[i]);
	}
	double*Vi = (double*)malloc(max*sizeof(double));
	for(i = 0;i < max;i++){
		Vi[i] = 0;
	}
	getInitial(Vi,C,H,O,N,type);
	//puts("0");
	double*Vii = newVector(Vi,P,K,C,H,O,N,size,type,I);
	/*for(i = 0;i < max;i++){
		printf("%f\n",Vii[i]);
	}printf("\n");getchar();*/
	/**///exit(1);
	for(i = 0;i < 7;i++){
		//printf("%.11f\n",K[i]);
	}//exit(1);
	int c = 0;
	while(/*errorV(Vii,Vi,size,Error,type)8*/c < lim){
		/*for(i = 0;i < max;i++){
		printf("%f\n",Vii[i]);
	}printf("\n----------------------\n");getchar();*/
		Vi = Vii;
		Vii = newVector(Vi,P,K,C,H,O,N,size,type,I);
		/*printf("%d\n",c+1);
		for(i = 0;i < max;i++){
		printf("%f\n",Vii[i]);
	}printf("\n");getchar();*/
		c++;
	}
	/*for(i = 0;i < max;i++){
		printf("%f\n",Vii[i]);
	}printf("\n");*///exit(1);
	completeProducts(p,Vii);
	
}
