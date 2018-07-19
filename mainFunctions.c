#include "mainFunctions.h"

bool Realloc(CReLIST*p){
	int i;
	CHEMICAL*aux;
	p->maxSize *= 2;
	aux = (CHEMICAL*)malloc(p->maxSize*sizeof(CHEMICAL));
	if(!aux){
		return True;
	}
	for(i = 0;i < p->count;i++){
		aux[i] = (p->List)[i];
	}
	free(p->List);
	p->List = aux;
	return False;
}

void ReallocCharMatrix(MATRIX*matrix,int size){
	int i;
	char**aux;
	aux = (char**)malloc(size*sizeof(char*));
	if(!aux){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	for(i = 0;i < matrix->count;i++){
		aux[i] = matrix->matrix[i];
	}
	free(matrix->matrix);
	matrix->matrix = aux;
}

void breakChemical(MATRIX*matrix,char component[]){
	int size = 20;
	int i,j,before,after,countInner = 0;
	int flag = 1;
	matrix->matrix = (char**)malloc(size*sizeof(char*));
	if(!matrix->matrix){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	matrix->count = 0;
	before = 0;
	after = 0;
	for(i=0;component[i];i++){
		if('0'<=component[i] && component[i]<='9' && flag){
			after = i;
			matrix->matrix[matrix->count] = (char*)malloc((after-before+1)*sizeof(char));
			if(!matrix->matrix[matrix->count]){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			for(j = before;j < after;j++){
				matrix->matrix[matrix->count][countInner] = component[j];
				countInner++;
			}
			matrix->matrix[matrix->count][countInner] = '\0';
			matrix->count++;
			countInner = 0;
			if(matrix->count==size){
				size*=2;
				ReallocCharMatrix(matrix,size);
			}
			before = after;
			flag = 0;
		}
		if('A'<=component[i] && component[i]<='Z' && !flag){
			after = i;
			matrix->matrix[matrix->count] = (char*)malloc((after-before+1)*sizeof(char));
			if(!matrix->matrix[matrix->count]){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			for(j = before;j<after;j++){
				matrix->matrix[matrix->count][countInner] = component[j];
				countInner++;
			}
			matrix->matrix[matrix->count][countInner] = '\0';
			matrix->count++;
			countInner = 0;
			if(matrix->count == size){
				size *= 2;
				ReallocCharMatrix(matrix,size);
			}
			before = after;
			flag = 1;
		}
	}
	after = i;
	matrix->matrix[matrix->count] = (char*)malloc((after-before+1)*sizeof(char));
	if(!matrix->matrix[matrix->count]){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	for(j = before;j<after;j++){
		matrix->matrix[matrix->count][countInner] = component[j];
		countInner++;
	}
	matrix->matrix[matrix->count][countInner] = '\0';
	matrix->count++;
}

ATOM fillComponents(MATRIX matrix){
	int len = matrix.count;
	int i,j;
	float value;
	char lastValue;
	ATOM p;
	p.C = 0;
	p.H = 0;
	p.O = 0;
	p.N = 0;
	p.S = 0;
	for(i = 0;i < len;i++){
		if(i%2){
			sscanf(matrix.matrix[i],"%f",&value);
			switch(lastValue){
				case 'C':
					p.C += value - 1;
					break;
				case 'H':
					p.H += value - 1;
					break;
				case 'O':
					p.O += value - 1;
					break;
				case 'N':
					p.N += value - 1;
					break;
				case 'S':
					p.S += value - 1;
					break;
				default:
					puts("Invalid atom in composition.");
					exit(1);
			}
		}
		else{
			for(j=0;j<strlen(matrix.matrix[i]);j++){
				lastValue = matrix.matrix[i][j];
				switch(lastValue){
					case 'C':
						p.C += 1;
						break;
					case 'H':
						p.H += 1;
						break;
					case 'O':
						p.O += 1;
						break;
					case 'N':
						p.N += 1;
						break;
					case 'S':
						p.S += 1;
						break;
					default:
						puts("Invalid atom in composition.");
						exit(1);
				}
			}
		}
	}
	return p;
}			

bool Collect(CReLIST*p,int flag){
	double count = 0;
	MATRIX matrix;
	char component[1000];
	
	while(count < 0.99995 || flag){
		if(!flag){
			if(!p->count){
				printf("\nInsert Fuel:\n\n");
			}
			else if(p->count == 1){
				printf("\nInsert %dnd Fuel:\n\n",p->count+1);
			}
			else if(p->count == 2){
				printf("\nInsert %drd Fuel:\n\n",p->count+1);
			}
			else{
				printf("\nInsert %dth Fuel:\n\n",p->count+1);
			}
		}
		else{
			if(!p->count){
				printf("\nInsert Product:\n\n");
			}
			else if(p->count == 1){
				printf("\nInsert %dnd Product:\n\n",p->count+1);
			}
			else if(p->count == 2){
				printf("\nInsert %drd Product:\n\n",p->count+1);
			}
			else{
				printf("\nInsert %dth Product:\n\n",p->count+1);
			}
		}
		p->List[p->count].got = 0;
		scanf("%lf",&p->List[p->count].Mol);
		if(p->List[p->count].Mol < 0 && flag){
			break;
		}
		getchar();
		fgets(component,1000,stdin);
		component[strlen(component) - 1] = '\0';
		count += p->List[p->count].Mol;
		
		p->List[p->count].Name = (char*)malloc(strlen(component)*sizeof(char));
		if(!p->List[p->count].Name){
			puts("Memory allocation error.\nTerminating the program.\n");
			exit(1);
		}
		
		strcpy(p->List[p->count].Name,component);

		breakChemical(&matrix,component);
		p->List[p->count].atoms = fillComponents(matrix);
		p->count++;
		
		if(p->count == p->maxSize){
			if(Realloc(p)){
				return True;
			}
		}
		
		if(count>1 && !flag){
			printf("Number of Mols in mixture surpassed the limit of 1 Mol. Restart.\n\n");
			count = 0;
			p->count = 0;
		}
	}
	free(matrix.matrix);
	return False;
}

double generateTeoricalAir(CReLIST r){
	int i;
	double NiC = 0;
	double NiO = 0;
	double NiH = 0;
	for(i = 0;i < r.count;i++){
		NiC += r.List[i].Mol*r.List[i].atoms.C;
		NiH += r.List[i].Mol*r.List[i].atoms.H;
		NiO -= r.List[i].Mol*r.List[i].atoms.O;
	}
	NiH /= 2;
	NiO += (2*NiC + NiH);
	NiO /= 2;
	
	return NiO;
}

void Combustion(CReLIST r,CReLIST*p,char c,double*air,double*nFuel,double*C,int P){
	int i;
	double NiCcomplete = 0, NiHcomplete = 0, NiNcomplete = 0,NiScomplete = 0,NiOcomplete = 0;
	double NiCincomplete = 0, NiHincomplete = 0,NiOincomplete = 0;
	double NiFuel,NiH2O,NiN2 = 100;
	for(i = 0;i < r.count;i++){
		NiCcomplete += r.List[i].Mol*r.List[i].atoms.C;
		NiHcomplete += r.List[i].Mol*r.List[i].atoms.H;
		NiNcomplete += r.List[i].Mol*r.List[i].atoms.N;
		NiScomplete += r.List[i].Mol*r.List[i].atoms.S;
		NiOcomplete += r.List[i].Mol*r.List[i].atoms.O;
	}
	*C = NiCcomplete;
	switch(c){
		case 'Y':
		case 'y':
		
			p->List[0].atoms.C = 1;
			p->List[0].atoms.H = 0;
			p->List[0].atoms.O = 2;
			p->List[0].atoms.N = 0;
			p->List[0].atoms.S = 0;
			p->List[0].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[0].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[0].Name,"CO2");
			p->List[0].Name[3] = '\0';
			p->List[0].Mol = NiCcomplete;
			p->List[0].got = 0;
			
			p->List[1].atoms.C = 0;
			p->List[1].atoms.H = 2;
			p->List[1].atoms.O = 1;
			p->List[1].atoms.N = 0;
			p->List[1].atoms.S = 0;
			p->List[1].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[1].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[1].Name,"H2O");
			p->List[1].Name[3] = '\0';
			p->List[1].Mol = NiHcomplete/2;
			p->List[1].got = 0;
		
			p->List[2].atoms.C = 0;
			p->List[2].atoms.H = 2;
			p->List[2].atoms.O = 0;
			p->List[2].atoms.N = 0;
			p->List[2].atoms.S = 0;
			p->List[2].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[2].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[2].Name,"H2");
			p->List[2].Name[3] = '\0';
			p->List[2].Mol = 0;
			p->List[2].got = 0;
		
			p->List[3].atoms.C = 0;
			p->List[3].atoms.H = 1;
			p->List[3].atoms.O = 0;
			p->List[3].atoms.N = 0;
			p->List[3].atoms.S = 0;
			p->List[3].Name = (char*)malloc(2*sizeof(char));
			if(!p->List[3].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[3].Name,"H");
			p->List[3].Name[1] = '\0';
			p->List[3].Mol = 0;
			p->List[3].got = 0;
		
			p->List[4].atoms.C = 0;
			p->List[4].atoms.H = 0;
			p->List[4].atoms.O = 2;
			p->List[4].atoms.N = 0;
			p->List[4].atoms.S = 0;
			p->List[4].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[4].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[4].Name,"O2");
			p->List[4].Name[2] = '\0';
			p->List[4].Mol = 0;
			p->List[4].got = 0;
		
			p->List[5].atoms.C = 0;
			p->List[5].atoms.H = 0;
			p->List[5].atoms.O = 0;
			p->List[5].atoms.N = 2;
			p->List[5].atoms.S = 0;
			p->List[5].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[5].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[5].Name,"N2");
			p->List[5].Name[2] = '\0';
			p->List[5].Mol = (NiNcomplete + P*7.52*(*air))/2;
			p->List[5].got = 0;
			
			p->List[6].atoms.C = 0;
			p->List[6].atoms.H = 0;
			p->List[6].atoms.O = 2;
			p->List[6].atoms.N = 1;
			p->List[6].atoms.S = 0;
			p->List[6].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[6].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[6].Name,"NO2");
			p->List[6].Name[3] = '\0';
			p->List[6].Mol = 0;
			p->List[6].got = 0;
			
			p->List[7].atoms.C = 0;
			p->List[7].atoms.H = 0;
			p->List[7].atoms.O = 2;
			p->List[7].atoms.N = 0;
			p->List[7].atoms.S = 1;
			p->List[7].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[7].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[7].Name,"SO2");
			p->List[7].Name[3] = '\0';
			p->List[7].Mol = NiScomplete;
			p->List[7].got = 0;
			
			p->List[8].atoms.C = 1;
			p->List[8].atoms.H = 0;
			p->List[8].atoms.O = 1;
			p->List[8].atoms.N = 0;
			p->List[8].atoms.S = 0;
			p->List[8].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[8].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[8].Name,"CO");
			p->List[8].Name[2] = '\0';
			p->List[8].Mol = 0;
			p->List[8].got = 0;
			
			p->List[9].atoms.C = 0;
			p->List[9].atoms.H = 0;
			p->List[9].atoms.O = 1;
			p->List[9].atoms.N = 1;
			p->List[9].atoms.S = 0;
			p->List[9].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[9].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[9].Name,"NO");
			p->List[9].Name[2] = '\0';
			p->List[9].Mol = 0;
			p->List[9].got = 0;
			
			p->List[10].atoms.C = 0;
			p->List[10].atoms.H = 2;
			p->List[10].atoms.O = 1;
			p->List[10].atoms.N = 0;
			p->List[10].atoms.S = 0;
			p->List[10].Name = (char*)malloc(2*sizeof(char));
			if(!p->List[10].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[10].Name,"O");
			p->List[10].Name[1] = '\0';
			p->List[10].Mol = 0;
			p->List[10].got = 0;
			
			p->List[11].atoms.C = 0;
			p->List[11].atoms.H = 1;
			p->List[11].atoms.O = 1;
			p->List[11].atoms.N = 0;
			p->List[11].atoms.S = 0;
			p->List[11].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[11].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[11].Name,"OH");
			p->List[11].Name[2] = '\0';
			p->List[11].Mol = 0;
			p->List[11].got = 0;
			
			p->count = 12;
			break;
			
		case 'N':
		case 'n':
			Collect(p,1);
			for(i = 0; i < p->count;i++){
				NiCincomplete += p->List[i].Mol*p->List[i].atoms.C;
				NiHincomplete += p->List[i].Mol*p->List[i].atoms.H;
			}
			NiFuel = NiCincomplete/NiCcomplete;
			*nFuel = NiFuel;
			NiOcomplete *= NiFuel;
			NiH2O = (NiFuel*NiHcomplete - NiHincomplete)/2;
			for(i = 0;i < p->count;i++){
				NiN2 -= p->List[i].Mol;
			}
			if(p->count > p->maxSize - 3){
				Realloc(p);
			}
			p->List[p->count].atoms.C = 0;
			p->List[p->count].atoms.H = 2;
			p->List[p->count].atoms.O = 1;
			p->List[p->count].atoms.N = 0;
			p->List[p->count].atoms.S = 0;
			p->List[p->count].Mol = NiH2O;
			p->List[p->count].Name = (char*)malloc(4*sizeof(char));
			if(!p->List[p->count].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[p->count].Name,"H2O");
			p->List[p->count].Name[3] = '\0';
			p->count++;
			
			p->List[p->count].atoms.C = 0;
			p->List[p->count].atoms.H = 0;
			p->List[p->count].atoms.O = 0;
			p->List[p->count].atoms.N = 2;
			p->List[p->count].atoms.S = 0;
			p->List[p->count].Mol = NiN2;
			p->List[p->count].Name = (char*)malloc(3*sizeof(char));
			if(!p->List[p->count].Name){
				puts("Memory allocation error.\nTerminating the program.\n");
				exit(1);
			}
			strcpy(p->List[p->count].Name,"N2");
			p->List[p->count].Name[2] = '\0';
			p->count++;
			break;
	}
	for(i = 0;i < p->count;i++){
		NiOincomplete += p->List[i].Mol*p->List[i].atoms.O;
	}
	*air = (NiOincomplete - NiOcomplete)/2;
}

double formula(Cp cp,double T){
	return R*( -cp.a1/T + cp.a2*log(T) + cp.a3*T + (cp.a4*pow(T,2))/2 + (cp.a5*pow(T,3))/3 + (cp.a6*pow(T,4))/4 + (cp.a7*pow(T,5))/5 + cp.b1);
}

double generateDh(CHEMICAL c,int flag,double n){
	double dh;
	//puts(c.Name);
	if(flag){
		dh = n*c.Mol*formula(c.cpLowRange,c.Temperature);
	}
	else{
		dh = n*c.Mol*formula(c.cpHighRange,c.Temperature);
	}
	//printf("%f\n",dh);
	return dh;
}

int estimateTemperature(CReLIST Reagents,CReLIST Products,CHEMICAL O2Air,CHEMICAL N2Air,double tAir,double nFuel,double*deltaH,int flag){
	int i;
	double energyR = 0,energyP = 0;
	double dh;
	float mols = 0;
	
	for(i = 0;i < Reagents.count;i++){
		if(Reagents.List[i].Temperature >= 200 && Reagents.List[i].Temperature < 1000){
			dh = generateDh(Reagents.List[i],1,nFuel);
		}
		if(Reagents.List[i].Temperature >= 1000 && Reagents.List[i].Temperature <= 6000){
			dh = generateDh(Reagents.List[i],0,nFuel);
		}
		energyR += dh;
	}
	for(i = 0;i < Products.count;i++){
		mols += Products.List[i].Mol;
		energyP += Products.List[i].Hf0 * Products.List[i].Mol;
		//puts(Products.List[i].Name);
		//printf("dh: %f\n",Products.List[i].Hf0 * Products.List[i].Mol);
	}
	
	if(O2Air.Temperature >= 200 && O2Air.Temperature < 1000){
		dh = generateDh(O2Air,1,1);
		energyR += tAir*dh;
		dh = generateDh(N2Air,1,1);
		energyR += tAir*dh;
	}
	if(O2Air.Temperature >= 1000 && O2Air.Temperature <= 6000){
		dh = generateDh(O2Air,0,1);
		energyR += tAir*dh;
		dh = generateDh(N2Air,0,1);
		energyR += tAir*dh;
	}
	*deltaH = energyR;
	dh = (energyR - energyP)/mols;
	printf("dh: %f\n",dh);
	
	if(flag){
		int h1,h2;
		int t1,t2;
		char array[50];
		FILE*arq;
	
		arq = fopen("table2.txt","r");
		fgets(array,50,arq);
		sscanf(array,"%d %d",&t1,&h1);
		fgets(array,50,arq);
		sscanf(array,"%d %d",&t2,&h2);
		while(1){
			if(dh<=h2){
				break;
			}
			if(t2 < 0){
				break;
			}
			t1 = t2; h1 = h2;
			fgets(array,50,arq);
			sscanf(array,"%d %d",&t2,&h2);
		}
		fclose(arq);
		return t1;
	}
	return 0;
}
	
void getValues(CReLIST*r,CReLIST*p){
	FILE*arq;
	char array[100];
	char name[20];
	int i;
	arq = fopen("table3.txt","r");
	while(fgets(array,100,arq)){
		sscanf(array,"%s",name);
		for(i = 0;i < r->count;i++){
			if(!strcmp(r->List[i].Name,name) && !(r->List[i].got)){
				r->List[i].got = 1;
				fgets(array,100,arq);
				fgets(array,100,arq);
				fgets(array,100,arq);
				sscanf(array,"%lf %lf %lf %lf %lf",&r->List[i].cpLowRange.a1,&r->List[i].cpLowRange.a2,&r->List[i].cpLowRange.a3,&r->List[i].cpLowRange.a4,&r->List[i].cpLowRange.a5);
				fgets(array,100,arq);
				sscanf(array,"%lf %lf %lf",&r->List[i].cpLowRange.a6,&r->List[i].cpLowRange.a7,&r->List[i].cpLowRange.b1);
				fgets(array,100,arq);
				fgets(array,100,arq);
				sscanf(array,"%lf %lf %lf %lf %lf",&r->List[i].cpHighRange.a1,&r->List[i].cpHighRange.a2,&r->List[i].cpHighRange.a3,&r->List[i].cpHighRange.a4,&r->List[i].cpHighRange.a5);
				fgets(array,100,arq);
				sscanf(array,"%lf %lf %lf",&r->List[i].cpHighRange.a6,&r->List[i].cpHighRange.a7,&r->List[i].cpHighRange.b1);
				rewind(arq);
			}
		}
	}
	rewind(arq);
	while(fgets(array,100,arq)){
		sscanf(array,"%s",name);
		for(i = 0;i < p->count;i++){
			if(!strcmp(p->List[i].Name,name) && !(p->List[i].got)){
				//puts(array);
				//puts(p->List[i].Name);
				p->List[i].got = 1;
				fgets(array,100,arq);
				fgets(array,100,arq);
				fgets(array,100,arq);
				//puts(array);
				sscanf(array,"%lf %lf %lf %lf %lf",&p->List[i].cpLowRange.a1,&p->List[i].cpLowRange.a2,&p->List[i].cpLowRange.a3,&p->List[i].cpLowRange.a4,&p->List[i].cpLowRange.a5);
				fgets(array,100,arq);
				//puts(array);
				sscanf(array,"%lf %lf %lf",&p->List[i].cpLowRange.a6,&p->List[i].cpLowRange.a7,&p->List[i].cpLowRange.b1);
				fgets(array,100,arq);
				fgets(array,100,arq);
				//puts(array);
				sscanf(array,"%lf %lf %lf %lf %lf",&p->List[i].cpHighRange.a1,&p->List[i].cpHighRange.a2,&p->List[i].cpHighRange.a3,&p->List[i].cpHighRange.a4,&p->List[i].cpHighRange.a5);
				fgets(array,100,arq);
				//puts(array);
				sscanf(array,"%lf %lf %lf",&p->List[i].cpHighRange.a6,&p->List[i].cpHighRange.a7,&p->List[i].cpHighRange.b1);
				rewind(arq);
				//printf("%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",p->List[i].cpHighRange.a1,p->List[i].cpHighRange.a2,p->List[i].cpHighRange.a3,p->List[i].cpHighRange.a4,p->List[i].cpHighRange.a5,p->List[i].cpHighRange.a6,p->List[i].cpHighRange.a7,p->List[i].cpHighRange.b1);
			}
		}
	}
	fclose(arq);
}	

void getHf0(CReLIST*p){
	FILE*arq;
	char array[100];
	char name[20];
	int i;
	double value;
	arq = fopen("table1.txt","r");
	while(fgets(array,100,arq)){
		sscanf(array,"%s %lf",name,&value);
		for(i = 0;i < p->count;i++){
			if(!strcmp(p->List[i].Name,name)){
				if(!(p->List[i].got)){
					p->List[i].Hf0 = 1000*value;
				}
			}
		}
	}
	fclose(arq);
}

double formulaIterativeD(double*v,double T){
	return R*((v[0]/pow(T,2)) + v[1]/T + v[2] + v[3]*T + v[4]*pow(T,2) + v[5]*pow(T,3) + v[6]*pow(T,4));
}

double formulaIterative(double*v,double T,double sum){
	return R*((-v[0]/T) + v[1]*log(T) + v[2]*T + v[3]*pow(T,2)/2 + v[4]*pow(T,3)/3 + v[5]*pow(T,4)/4 + v[6]*pow(T,5)/5 + v[7]) - sum;
}
void printaV(double*v,double sum){
	int i;
	printf("8.314*(\n");
	for(i = 0;i < 8;i++){
		if(!i){
			printf("%.20f/T\n",-v[i]);
		}
		else if(i == 1){
			printf("%.20f*log(T)\n",v[i]);
		}
		else if(i == 7){
			printf("%.20f\n",v[i]);
		}
		else{
			printf("%.20f*(T^%d)/%d\n",v[i],i - 1,i - 1);
		}
	}
	printf(") - %f\n",sum);
}

double newtonFormula(double*v,double T, double sum){
	//printf("teste: %f\n",sum);
	return T - formulaIterative(v,T,sum)/formulaIterativeD(v,T);
}

void calculateCp(double*v,double temp,CReLIST p){
	int i,j;
	for(j = 0;j < 8;j++){
		v[j] = 0;
	}
	for(i = 0;i < p.count;i++){
		//puts(p.List[i].Name);
		//printf("Mol: %f\n",p.List[i].Mol);
		if(200 <= temp && temp < 1000){
			for(j = 0;j < 8;j++){
				switch(j){
					case 0:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a1;
						break;
					case 1:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a2;
						break;
					case 2:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a3;
						break;
					case 3:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a4;
						break;
					case 4:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a5;
						break;
					case 5:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a6;
						break;
					case 6:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.a7;
						break;
					case 7:
						v[j] += p.List[i].Mol*p.List[i].cpLowRange.b1;
						break;
				}
			}
		}
		else{
			for(j = 0;j < 8;j++){
				switch(j){
					case 0:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a1;
						break;
					case 1:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a2;
						break;
					case 2:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a3;
						break;
					case 3:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a4;
						break;
					case 4:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a5;
						break;
					case 5:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a6;
						break;
					case 6:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.a7;
						break;
					case 7:
						v[j] += p.List[i].Mol*p.List[i].cpHighRange.b1;
						break;
				}
			}
		}
		int k;
	/*for(k = 0;k < 8;k++){
		printf("V[%d]: %.11f\n",k,v[k]);
	}
	printf("----------\n");*/
	}
}

double newtonIterativeMethod(double temp,double*v,double sum,CReLIST p){
	//printaV(v,sum);
	double Ti,Tii;
	int isHigh;
	int check = 0;
	if(temp >= 1000){
		isHigh = 1;
	}
	else{
		isHigh = 0;
	}
	Ti = /*(double)*/temp;
	//printf("Antes: %f\n",Ti);
	Tii = newtonFormula(v,Ti,sum);
	//Tii = newtonFormula(v,temp,sum);
	//Tii = Bissection(v,temp,sum);
	//printf("Depois: %f\n",Tii);
	while(abs(Tii - Ti) >= 0.1){
		//printf("i: %d\n",i);i++;
		if((isHigh && Tii<1000)||(!isHigh && Tii>=1000)){
			calculateCp(v,Tii,p);
			isHigh = !isHigh;
		}
		Ti = Tii;
		Tii = newtonFormula(v,Ti,sum);
		check++;
		if(check == 1000){
			puts("Iterative method reached endless loop.\nTerminating the application.\n");
			exit(1);
		}
	}
	return Tii;
}

double generateResults(double temp, CReLIST p,double dh){
	double sum = dh;
	double*vectorSumCp;
	double*vectorSumCp2;
	double result;
	vectorSumCp = (double*)malloc(8*sizeof(double));
	vectorSumCp2 = (double*)malloc(8*sizeof(double));
	if(!vectorSumCp){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}

	calculateCp(vectorSumCp,temp,p);calculateCp(vectorSumCp2,temp,p);
	result = newtonIterativeMethod(temp,vectorSumCp,sum,p);
	free(vectorSumCp);
	return result;
}

int saveData(char*name,CReLIST r, CReLIST p,CHEMICAL O2Air,double T,double n,double air){
	int i;
	time_t t = time(NULL);
	struct tm tm = *localtime(&t);
	FILE*arq;
	arq = fopen(name,"a");
	if(!arq){
		return 0;
	}
	fprintf(arq,"-------------------------------------------------------------\n");
	fprintf(arq,"						RESULTS (%d/%d/%d)\n",tm.tm_mday, tm.tm_mon + 1,tm.tm_year + 1900);
	fprintf(arq,"-------------------------------------------------------------\n");
	for(i = 0;i < r.count;i++){
		fprintf(arq,"%.6f %s + ",n*r.List[i].Mol,r.List[i].Name);
	}
	fprintf(arq,"%.6f(O2 + 3.76N2) -> ",air);
	for(i = 0;i < p.count;i++){
		if(i == p.count - 1){
			fprintf(arq,"%.6f %s\n\n",p.List[i].Mol,p.List[i].Name);
		}
		else{
			fprintf(arq,"%.6f %s + ",p.List[i].Mol,p.List[i].Name);
		}
	}
	for(i = 0; i < r.count;i++){
		fprintf(arq,"%s temperature: %fK\n",r.List[i].Name,r.List[i].Temperature);
	}
	fprintf(arq,"Air temperature: %fK\n",O2Air.Temperature);
	fprintf(arq,"\nTemperature estimative: %.6fK\n\n",T);
	fclose(arq);
	return 1;
}

void unGet(CReLIST*r,CReLIST*p){
	int i;
	for(i = 0;i < r->count;i++){
		r->List[i].got = 0;
	}
	for(i = 0;i < p->count;i++){
		p->List[i].got = 0;
	}
}

void resetProducts(CReLIST*p,double C,double H,double N){
	int i;
	for(i = 0;i < p->maxSize;i++){
		if(i != 7){
			p->List[i].Mol = 0;
		}
	}
	p->List[0].Mol = C;
	p->List[1].Mol = 2;
	p->List[5].Mol = N;
	p->List[2].Mol = 2;
}

