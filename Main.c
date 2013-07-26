#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "InitiationFunctions.h"
#include "MutationFunctions.h"
#include "SelectionFunctions.h"
#include "MeiosisFunctions.h"
#include "SyngamyFunctions.h"
#include "BiogenesisFunctions.h"
#include "EnvironmentalFluctuations.h"
#include "DevelopmentFunctions.h"

void WriteGenotypeFile( double    mat[3][M+1], int choice, int index){
	char name[15];
	int n;
	
	FILE *pfile=NULL;
	
	switch (choice) {
		case 1:
			n=sprintf(name, "a_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 2:
			n=sprintf(name, "b_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
			
		default:
			n=sprintf(name, "c_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
	}
	
	int col;
	for (col=0; col<(M+1); col++)
		fprintf(pfile, "%.10f	%.10f	%.10f\n", mat[0][col], mat[1][col], mat[2][col]);
	fclose(pfile);
}

void WriteMemoryFile( double    mat[runs][3], int choice, int index){
	char name[15];
	int n;
	
	FILE *pfile=NULL;
	switch (choice) {
		case 1:
			n=sprintf(name, "memory_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 2:
			n=sprintf(name, "fitness_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
			
		default:
			n=sprintf(name, "variance_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
	}

	int row;
	for (row=0; row<(runs+1); row++)
		fprintf(pfile, "%.20f	%.20f	%.20f\n", mat[row][0], mat[row][1], mat[row][2]);
	fclose(pfile);
}


void WriteHomFile( double  mat[runs][1], int index){
	char name[15];
	int n;
	
	FILE *pfile=NULL;
	
	n=sprintf(name, "hom_%d.txt", index);
	pfile = fopen(name, "w+");
	
	int row;
	for (row=0; row<(runs+1); row++)
		fprintf(pfile, "%.10f\n", mat[row][0]);
	fclose(pfile);
}
		

//////////////////////////////////////////////////////////////////////////////
////////////////////////// M A I N   P R O G R A M //////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void) {
	/////////////////////////////////////
	//initialize matrices for each step// 
	////////////////////////////////////
	double  pm1[M+1][M+1];		//mutation matrix
	InitMutMat(pm1);
	double  bioprob[M+1][M+1];	//biogenesis probabilities
	InitBioMat(bioprob);
	double  pmeio1[M+1][M+1];	//first meiotic step
	InitMeio1(pmeio1);
	double  pmeio2[M+1][M+1];	//second meiotic step 
	InitMeio2(pmeio2);
	double    aa[NUC][M+1];		//mitonuclear state for aa
	InitMatrix(aa);
	InitMatrixa(aa);
	double    AA[NUC][M+1];		//mitonuclear state for AA
	InitMatrix(AA);  
//	InitMatrixa(AA);
	double    Aa[NUC][M+1];		//mitonuclear state for Aa 
	InitMatrix(Aa); 
	double    aagam[2][M+1];		// aa's gametes
	InitGameteMatrix(aagam);
	double    AAgam[2][M+1];		// AA's gametes
	InitGameteMatrix(AAgam);
	double    Aagam[2][M+1];		// Aa's gametes
	InitGameteMatrix(Aagam);
	double    agamete[2][M+1];		// gametes with gene a
	InitGameteMatrix(agamete);
	double    Agamete[2][M+1];		// gametes with gene A
	InitGameteMatrix(Agamete);
	double    faa[NUC][M+1];		//non-changing matrix used for aa
	InitMatrix(faa);   	
	double    fAA[NUC][M+1];		//non-changing matrix used for AA
	InitMatrix(fAA);
	double    fAa[NUC][M+1];		//non-changing matrix used for Aa
	InitMatrix(fAa);
	double    memory[runs][3];		//matrix for storing totals in each group 
	InitMem(memory);
	double    fitnessmem[runs][3];	//matrix for storing average in each group 
	InitMem(fitnessmem);
	double    variancemem[runs][3];	//matrix for storing variance in each group 
	InitMem(variancemem);
	double    hom[runs][1];		//matric for storing the dominant homozygote state
	InitHom(hom);
	
	//print initial totals for a,b and c
	printf("\n");
	printf("%.10f0", SumMatrix(aa));
	printf("\n");
	printf("%.10f0", SumMatrix(AA));
	printf("\n");
	printf("%.10f0", SumMatrix(Aa));
	printf("\n");	
	printf("\n");	
	
	
	/////////////////////////////////////
	////// Main loop - Life cycle/////// 
	///////////////////////////////////
	int index;
	double    homozygote=0.0;
	
	for (index=0; index<(runs); index++) {

		//hom[index][0] = SumNuc(aa, 0) + SumNuc(Aa, 0) + SumNuc(AA, 0);
		//if (index%changenuc==0) {
		//	homozygote = hom[index][0];
		//}
		
		//introduce mutation at step 100
		if (index==100) {
			InitMatrixc(aa, Aa, 100); 
			
		}
		
		if (index==102 || index==2000) {
			printf("\n");
			printf("%.10f0", SumMatrix(aa));
			printf("\n");
			printf("%.10f0", SumMatrix(AA));
			printf("\n");
			printf("%.10f0", SumMatrix(Aa));
			printf("\n");	
			printf("\n");		
			printf("%.10f0", memory[index-1][2]/2+ memory[index-1][1]);
			printf("\n");
		}
	
		
		memory[index][0] = SumMatrix(aa);
		memory[index][1] = SumMatrix(AA);
		memory[index][2] = SumMatrix(Aa); 
		
		fitnessmem[index][0] = wbarType(aa);
		fitnessmem[index][1] = wbarType(AA);
		fitnessmem[index][2] = wbarType(Aa);
		
		variancemem[index][0] = VarType(aa);
		variancemem[index][1] = VarType(AA);
		variancemem[index][2] = VarType(Aa);
		
		int col;
		/////////////
		//MUTATION// 
		///////////
		CopyMat(aa,faa);
		CopyMat(AA,fAA);
		CopyMat(Aa,fAa);
		
		for (col=0; col<(M+1); col++) {
			Mutation(pm1, faa, aa, 0, col);
			Mutation(pm1, faa, aa, 1, col);
			Mutation(pm1, faa, aa, 2, col);
			
			Mutation(pm1, fAA, AA, 0, col);
			Mutation(pm1, fAA, AA, 1, col);
			Mutation(pm1, fAA, AA, 2, col);
			
			Mutation(pm1, fAa, Aa, 0, col);
			Mutation(pm1, fAa, Aa, 1, col);
			Mutation(pm1, fAa, Aa, 2, col);
		}
		
		/////////////////////////////		
		//MITOCHONDRIAL BIOGENESIS//
		///////////////////////////
		//CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
				
		//Biogenesis(faa, aa, bioprob);
		
		//Biogenesis(fAA, AA, bioprob);
		
		//Biogenesis(fAa, Aa, bioprob);
		
		//////////////
		//SELECTION//
		/////////////
		//if (index<300) {
		//	Select(aa, AA, Aa); //only with a and A
		//}
		//else {
		//	SelectEnv(faa, aa, fAA,AA , fAa, Aa, index, homozygote);
		//}

		Select(aa, AA, Aa); //only with a and A
		//SelectEnv(faa, aa, fAA,AA , fAa, Aa, index, homozygote);

		//////////////		
		//Meiosis I//
		////////////
		CopyMat(aa,faa);
		CopyMat(AA,fAA);
		CopyMat(Aa,fAa);
		for (col=0; col<(M+1); col++) {
			Meio1(pmeio1, faa, aa, 0, col);
			Meio1(pmeio1, faa, aa, 1, col);
			Meio1(pmeio1, faa, aa, 2, col);
			
			Meio1(pmeio1, fAA, AA, 0, col);
			Meio1(pmeio1, fAA, AA, 1, col);
			Meio1(pmeio1, fAA, AA, 2, col);
			
			Meio1(pmeio1, fAa, Aa, 0, col);
			Meio1(pmeio1, fAa, Aa, 1, col);
			Meio1(pmeio1, fAa, Aa, 2, col);
		}

		///////////////
		//Meiosis II//
		/////////////
		CopyMat(aa,faa);
		CopyMat(AA,fAA);
		CopyMat(Aa,fAa);
		for (col=0; col<(M/2+1); col++) {
			Meio2(pmeio2, faa, aagam, 0, col);
			Meio2(pmeio2, faa, aagam, 1, col);
			
			Meio2(pmeio2, fAA, AAgam, 0, col);
			Meio2(pmeio2, fAA, AAgam, 1, col);
			
			Meio2(pmeio2, fAa, Aagam, 0, col);
			Meio2(pmeio2, fAa, Aagam, 1, col);
		}	

		//Define each gamete matrix
		GetGametes(aagam, AAgam, Aagam, agamete, Agamete); //when only a and A 
				
		//Syngamy
		for (col=0; col<(M+1); col++) {
			syngamy2(0, col, agamete, agamete, aa);
			syngamy2(1, col, agamete, agamete, aa);
			syngamy2(2, col, agamete, agamete, aa);
			
//			syngamy2(0, col, Agamete, Agamete, AA);
//			syngamy2(1, col, Agamete, Agamete, AA);
//			syngamy2(2, col, Agamete, Agamete, AA);
			
			syngamy1AxA(0, col, Agamete, Agamete, AA);
			syngamy1AxA(1, col, Agamete, Agamete, AA);
			syngamy1AxA(2, col, Agamete, Agamete, AA);
			
			syngamy1(0, col, Agamete, agamete, Aa);
			syngamy1(1, col, Agamete, agamete, Aa);
			syngamy1(2, col, Agamete, agamete, Aa);
		}		

	}
	//write data into file
//	WriteMemoryFile(memory, 1, index);
//	WriteMemoryFile(fitnessmem, 2, index);
//	WriteMemoryFile(variancemem, 3, index);
//	WriteGenotypeFile(aa, 1, index);
//	WriteGenotypeFile(AA, 2, index);
//	WriteGenotypeFile(Aa, 3, index);
//	WriteHomFile(hom, index);
	

	
	return 0;
	
}



