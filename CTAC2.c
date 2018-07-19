#include "init.h"
#include "mainFunctions.h"
#include "globals.h"
int main(void){
	
	char*fileName;
	
	char array[100];
	
	double Ti;
	
	double Tii;
	
	CReLIST Reagents, Products;
	
	CHEMICAL O2Air, N2Air;
	
	double deltaH;
	
	double C; double H;double O; double N;
	
	double P;
	
	double teoricalAir;
	
	double nFuel = 1;
	
	char isComplete = '\0';
	
	int pureO2;
	
	int i;
	printf("\n\n\n\n\n");
	
	printf("________________________________________________________________________________\n|                                                                              ||                    ADIABATIC FLAME TEMPERATURE CALCULATOR                    ||                                                                              |\n________________________________________________________________________________\n");
	
	printf("Burn with Theorical Air or O2 [1/0]\n\n");
	scanf("%d",&pureO2);
	
	getchar();
	while(pureO2 < 0 &&pureO2 > 1){
		scanf("%d",&pureO2);
		getchar();
	}
	
	initAir(&O2Air,&N2Air,pureO2);
	
	init_Matrix_of_Functions();

	if(init(&Reagents)){
		printf("Memory Error\n");
		return 1;
	}

	if(init(&Products)){
		printf("Memory Error\n");
		return 1;
	}
	
	if(Collect(&Reagents,0)){
		printf("Memory Error\n");
		return 1;
	}
	
	
	teoricalAir = generateTeoricalAir(Reagents);
	printf("\nPressure of the reaction in atm\n\n");
	scanf("%lf",&P);
	getchar();
	while(P < 0){
		puts("Invalid Pressure. Pressure must be a postive value\nPressure of the reaction in atm\n\n");
		scanf("%lf",&P);
		getchar();
	}
	printf("\nComplete Combustion Process?\n(Y/N)\n\n");
	isComplete = getchar();
	while(isComplete != 'Y' && isComplete != 'N' && isComplete != 'y' && isComplete != 'n'){
		printf("Invalid command!\nComplete Combustion Process?\n(Y/N)\n\n");
		isComplete = getchar();
	}
	
	Combustion(Reagents,&Products,isComplete,&teoricalAir,&nFuel,&C,pureO2);
	
	getValues(&Reagents,&Products);
	
	unGet(&Reagents,&Products);
	
	getHf0(&Products);
	int countPrinted = 0;
	printf("________________________________________________________________________________|                                                                              ||                          BALANCED CHEMICAL EQUATION                          ||                                                                              |________________________________________________________________________________\n");
	for(i = 0;i < Reagents.count;i++){
		printf("%.2f %s + ",nFuel*Reagents.List[i].Mol,Reagents.List[i].Name);
	}
	printf("%.2f(O2 + %.2fN2) -> ",teoricalAir,pureO2*3.76);
	for(i = 0;i < Products.count;i++){
		if(Products.List[i].Mol){
			if(!countPrinted){
				printf("%.2f %s",Products.List[i].Mol,Products.List[i].Name);
				countPrinted++;
			}
			else{
				printf(" + %.2f %s ",Products.List[i].Mol,Products.List[i].Name);
			}
		}
	}
	printf("\n________________________________________________________________________________");
	printf("\n\n");
	
	for(i = 0;i < Reagents.count;i++){
		printf("Enter Temperature for %s in Kelvin\n\n",Reagents.List[i].Name);
		scanf("%lf",&Reagents.List[i].Temperature);
		getchar();
	}
	
	printf("\nEnter Temperature for Air in Kelvin\n\n");
	scanf("%lf",&O2Air.Temperature);
	getchar();
	N2Air.Temperature = O2Air.Temperature;
	//Reagents.List[0].Mol = 4;O2Air.Mol = 1;
	Ti = estimateTemperature(Reagents,Products,O2Air,N2Air,teoricalAir,nFuel,&deltaH,1);
	printf("\nThe first estimative for the Products Temperature is: %fK\n\n",Ti);
	//-----------------------------
	Tii = generateResults(Ti,Products,deltaH);
	//printf("%f\n",Tii);exit(1);
	C = Products.List[0].Mol; H = Products.List[1].Mol; N = Products.List[5].Mol;O = 2*C + H + 2*teoricalAir;
	//---------------------------------------
	//printf("%f\n%f\n%f\n",C,H,N);exit(1);
	//int ccc = 0;
	H = 8;
	O = 2;
	while(abs(Tii - Ti) >= 0.01 && (isComplete == 'y' || isComplete == 'Y')){
		printf("T: %f\n",Tii);
		Ti = Tii;
		newProducts(&Products,Tii,P,C,H,O,N);
		Tii = generateResults(Ti,Products,deltaH);
	}
	newProducts(&Products,Tii,P,C,H,O,N);
	printf("The final estimative for the Products Temperature is: %.3fK\n\n",Tii);
	exit(1);
	printf("Specify the NAME of the FILE:\n\n");
	fgets(array,100,stdin);
	array[strlen(array) - 1] = '\0';
	fileName = (char*)malloc(strlen(array)*sizeof(char));
	strcpy(fileName,array);
	if(!fileName){
		puts("Memory allocation error.\nTerminating the program.\n");
		exit(1);
	}
	
	printf("\nSaving data in file\n");
	if(!saveData(fileName,Reagents,Products,O2Air,Tii,nFuel,teoricalAir)){
		puts("Error saving the data.\n");
	}
	else{
		puts("\nData saved successfully.\n");
	}
	
	free(Products.List);
	free(Reagents.List);
	
	printf("Press ENTER to return\n");
	getchar();
	
	return 0;
}
