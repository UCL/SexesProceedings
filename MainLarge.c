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


void WriteGenotypeFile( double  mat[3][M+1], int choice, int index){
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
		case 3:
			n=sprintf(name, "c_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 4:
			n=sprintf(name, "d_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 5:
			n=sprintf(name, "e_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 6:
			n=sprintf(name, "f_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		case 7:
			n=sprintf(name, "g_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
		default:
			n=sprintf(name, "h_%d.txt", index);
			pfile = fopen(name, "w+");
			break;
	}	
	int col;
	for (col=0; col<(M+1); col++)
		fprintf(pfile, "%.10f	%.10f	%.10f\n", mat[0][col], mat[1][col], mat[2][col]);
	fclose(pfile);
}

void WriteMemoryFile( double  mat[runs][8], int choice, int index){
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
		fprintf(pfile, "%.10f	%.10f	%.10f	%.10f	%.10f	%.10f	%.10f	%.10f\n", mat[row][0], mat[row][1], mat[row][2], mat[row][3], mat[row][4], mat[row][5], mat[row][6], mat[row][7]);
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

//Here we can introduce mutations to study the larger system with a, a1, a2, A1 and A2.
//Different matrices can be activated to study different scenarios. 

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
	 double  aa[NUC][M+1];		//mitonuclear state for aa
	InitMatrix(aa);
	InitMatrixa(aa);
	 double  AA[NUC][M+1];		//mitonuclear state for AA
	InitMatrix(AA); 
	 double  Aa[NUC][M+1];		//mitonuclear state for Aa 
	InitMatrix(Aa); 
	 double  A1a[NUC][M+1];		//mitonuclear state for A1a 
	InitMatrix(A1a);  
	 double  A1A[NUC][M+1];		//mitonuclear state for A1A 
	InitMatrix(A1A);  
	 double  a2a[NUC][M+1];		//mitonuclear state for a2a 
	InitMatrix(a2a);
//	InitMatrixa(a2a);
	 double  Aa2[NUC][M+1];		//mitonuclear state for Aa2 
	InitMatrix(Aa2);  
	 double  A1a2[NUC][M+1];		//mitonuclear state for A1a2 
	InitMatrix(A1a2);
	 double  aagam[2][M+1];		// aa's gametes
	InitGameteMatrix(aagam);
	 double  AAgam[2][M+1];		// AA's gametes
	InitGameteMatrix(AAgam);
	 double  Aagam[2][M+1];		// Aa's gametes
	InitGameteMatrix(Aagam);
	 double  A1agam[2][M+1];		// A1a's gametes
	InitGameteMatrix(A1agam);
	 double  A1Agam[2][M+1];		// A1A's gametes
	InitGameteMatrix(A1Agam);
	 double  a2agam[2][M+1];		// a2a's gametes
	InitGameteMatrix(a2agam);
	 double  Aa2gam[2][M+1];		// Aa2's gametes
	InitGameteMatrix(Aa2gam);
	 double  A1a2gam[2][M+1];		// A1a2's gametes
	InitGameteMatrix(A1a2gam);
	 double  agamete[2][M+1];		// gametes with gene a
	InitGameteMatrix(agamete);
	 double  Agamete[2][M+1];		// gametes with gene A
	InitGameteMatrix(Agamete);
	 double  a2gamete[2][M+1];		// gametes with gene a2
	InitGameteMatrix(a2gamete);
	 double  A1gamete[2][M+1];		// gametes with gene A1
	InitGameteMatrix(A1gamete);
	 double  faa[NUC][M+1];		//non-changing matrix used for aa
	InitMatrix(faa);   	
	 double  fAA[NUC][M+1];		//non-changing matrix used for AA
	InitMatrix(fAA);
	 double  fAa[NUC][M+1];		//non-changing matrix used for Aa
	InitMatrix(fAa);
	 double  fA1a[NUC][M+1];		//non-changing matrix used for A1a
	InitMatrix(fA1a);   	
	 double  fA1A[NUC][M+1];		//non-changing matrix used for A1A
	InitMatrix(fA1A);
	 double  fa2a[NUC][M+1];		//non-changing matrix used for a2a
	InitMatrix(fa2a);	
	 double  fAa2[NUC][M+1];		//non-changing matrix used for Aa2
	InitMatrix(fAa2);
	 double  fA1a2[NUC][M+1];		//non-changing matrix used for A1a2
	InitMatrix(fA1a2);	
	 double  memory[runs][8];		//matrix for storing totals in each group 
	InitMemLarge(memory);
	 double  fitnessmem[runs][8];	//matrix for storing average in each group 
	InitMemLarge(fitnessmem);
	 double  variancemem[runs][8];	//matrix for storing variance in each group 
	InitMemLarge(variancemem);
	 double  hom[runs][1];		//matric for storing the dominant homozygote state
	InitHom(hom);

	
	//print initial totals for a,b and c
	printf("\n");
	printf("%.10f0", SumMatrix(aa));
	printf("\n");
	printf("%.10f0", SumMatrix(AA));
	printf("\n");
	printf("%.10f0", SumMatrix(Aa));
	printf("\n");
	printf("A1a\n");
	printf("%.10f0", SumMatrix(A1a));
	printf("\n");
	printf("A1A\n");
	printf("%.10f0", SumMatrix(A1A));
	printf("\n");
	printf("a2a\n");
	printf("%.10f0", SumMatrix(a2a));
	printf("\n");
	printf("Aa2\n");
	printf("%.10f0", SumMatrix(Aa2));
	printf("\n");
	printf("A1a2\n");
	printf("%.10f0", SumMatrix(A1a2));
	printf("\n");
	printf("\n");
	
		
	/////////////////////////////////////
	////// Main loop - Life cycle/////// 
	///////////////////////////////////
	int index;
	 double  homozygote=0.0;
	
	for (index=0; index<runs; index++) {
		
		//keep track of homozygote state
		//hom[index][0] = SumNuc(aa, 0) + SumNuc(Aa, 0)+ SumNuc(AA, 0);
		//if (index%changenuc==0) {
		//	homozygote = hom[index][0];
		//}
		
		//introduce mutations
		if (index==100) InitMatrixc(aa, A1a, 100.0); 

		memory[index][0] = SumMatrix(aa);
		memory[index][1] = SumMatrix(AA);
		memory[index][2] = SumMatrix(Aa);
		memory[index][3] = SumMatrix(A1a);
		memory[index][4] = SumMatrix(A1A);
		memory[index][5] = SumMatrix(a2a);
		memory[index][6] = SumMatrix(Aa2);
		memory[index][7] = SumMatrix(A1a2);
		
		fitnessmem[index][0] = wbarType(aa);
		fitnessmem[index][1] = wbarType(AA);
		fitnessmem[index][2] = wbarType(Aa);
		fitnessmem[index][3] = wbarType(A1a);
		fitnessmem[index][4] = wbarType(A1A);
		fitnessmem[index][5] = wbarType(a2a);
		fitnessmem[index][6] = wbarType(Aa2);
		fitnessmem[index][7] = wbarType(A1a2);
		
		variancemem[index][0] = VarType(aa);
		variancemem[index][1] = VarType(AA);
		variancemem[index][2] = VarType(Aa);
		variancemem[index][3] = VarType(A1a);
		variancemem[index][4] = VarType(A1A);
		variancemem[index][5] = VarType(a2a);
		variancemem[index][6] = VarType(Aa2);
		variancemem[index][7] = VarType(A1a2);	
		
		int col;
		/////////////
		//MUTATION// 
		///////////
		CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
		CopyMat(A1a,fA1a);
		//CopyMat(A1A,fA1A);
		//CopyMat(a2a,fa2a);
		//CopyMat(Aa2,fAa2);
		CopyMat(A1a2,fA1a2);
		
		for (col=0; col<(M+1); col++) {
			Mutation(pm1, faa, aa, 0, col);
			Mutation(pm1, faa, aa, 1, col);
			Mutation(pm1, faa, aa, 2, col);
			
	//		Mutation(pm1, fAA, AA, 0, col);
	//		Mutation(pm1, fAA, AA, 1, col);
	//		Mutation(pm1, fAA, AA, 2, col);
			
	//		Mutation(pm1, fAa, Aa, 0, col);
	//		Mutation(pm1, fAa, Aa, 1, col);
	//		Mutation(pm1, fAa, Aa, 2, col);
			
			Mutation(pm1, fA1a, A1a, 0, col);
			Mutation(pm1, fA1a, A1a, 1, col);
			Mutation(pm1, fA1a, A1a, 2, col);
			
	//		Mutation(pm1, fA1A, A1A, 0, col);
	//		Mutation(pm1, fA1A, A1A, 1, col);
	//		Mutation(pm1, fA1A, A1A, 2, col);
			
	//		Mutation(pm1, fa2a, a2a, 0, col);
	//		Mutation(pm1, fa2a, a2a, 1, col);
	//		Mutation(pm1, fa2a, a2a, 2, col);
		
	//		Mutation(pm1, fAa2, Aa2, 0, col);
	//		Mutation(pm1, fAa2, Aa2, 1, col);
	//		Mutation(pm1, fAa2, Aa2, 2, col);
			
			Mutation(pm1, fA1a2, A1a2, 0, col);
			Mutation(pm1, fA1a2, A1a2, 1, col);
			Mutation(pm1, fA1a2, A1a2, 2, col);
		}
		
		/////////////////////////////		
		//MITOCHONDRIAL BIOGENESIS//
		///////////////////////////
		//need to check this works 
		//CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
		//CopyMat(A1a,fA1a);
		//CopyMat(A1A,fA1A);
		//CopyMat(a2a,fa2a);
		//CopyMat(Aa2,fAa2);
		//CopyMat(A1a2,fA1a2);
		
		//Biogenesis(faa, aa, bioprob);
				
		//Biogenesis(fAA, AA, bioprob);
		
		//Biogenesis(fAa, Aa, bioprob);
		
		//Biogenesis(fA1a, A1a, bioprob);
		
		//Biogenesis(fA1A, A1A, bioprob);
		
		//Biogenesis(fa2a, a2a, bioprob);
		
		//Biogenesis(fAa2, Aa2, bioprob);
		
		//Biogenesis(fA1a2, A1a2, bioprob);

		//////////////
		//SELECTION//
		/////////////
		CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
		CopyMat(A1a,fA1a);
		//CopyMat(A1A,fA1A);
		//CopyMat(a2a,fa2a);
		//CopyMat(Aa2,fAa2);
		CopyMat(A1a2,fA1a2);
		
//		if (index>2000) {
//			SelectionLarge(faa,aa,fAA, AA, fAa, Aa, fA1a, A1a, fA1A, A1A, fa2a, a2a, fAa2, Aa2, fA1a2, A1a2);
//		}
//		else {
//			SelectionLargeEnv(faa,aa,fAA, AA, fAa, Aa, fA1a, A1a, fA1A, A1A, fa2a, a2a, fAa2, Aa2, fA1a2, A1a2, index, homozygote);
//		}

		SelectionLarge(faa,aa,fAA, AA, fAa, Aa, fA1a, A1a, fA1A, A1A, fa2a, a2a, fAa2, Aa2, fA1a2, A1a2);
		//SelectionLargeEnv(faa,aa,fAA, AA, fAa, Aa, fA1a, A1a, fA1A, A1A, fa2a, a2a, fAa2, Aa2, fA1a2, A1a2, index, homozygote);

		//////////////		
		//Meiosis I//
		////////////
		CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
		CopyMat(A1a,fA1a);
		//CopyMat(A1A,fA1A);
		//CopyMat(a2a,fa2a);
		//CopyMat(Aa2,fAa2);
		CopyMat(A1a2,fA1a2);
		
		for (col=0; col<(M+1); col++) {
			Meio1(pmeio1, faa, aa, 0, col);
			Meio1(pmeio1, faa, aa, 1, col);
			Meio1(pmeio1, faa, aa, 2, col);
			
		//	Meio1(pmeio1, fAA, AA, 0, col);
		//	Meio1(pmeio1, fAA, AA, 1, col);
		//	Meio1(pmeio1, fAA, AA, 2, col);
			
		//	Meio1(pmeio1, fAa, Aa, 0, col);
		//	Meio1(pmeio1, fAa, Aa, 1, col);
		//	Meio1(pmeio1, fAa, Aa, 2, col);
			
			
			Meio1(pmeio1, fA1a, A1a, 0, col);
			Meio1(pmeio1, fA1a, A1a, 1, col);
			Meio1(pmeio1, fA1a, A1a, 2, col);
			
		//	Meio1(pmeio1, fA1A, A1A, 0, col);
		//	Meio1(pmeio1, fA1A, A1A, 1, col);
		//	Meio1(pmeio1, fA1A, A1A, 2, col);
			
		//	Meio1(pmeio1, fa2a, a2a, 0, col);
		//	Meio1(pmeio1, fa2a, a2a, 1, col);
		//	Meio1(pmeio1, fa2a, a2a, 2, col);
			
		//	Meio1(pmeio1, fAa2, Aa2, 0, col);
		//	Meio1(pmeio1, fAa2, Aa2, 1, col);
		//	Meio1(pmeio1, fAa2, Aa2, 2, col);
			
			Meio1(pmeio1, fA1a2, A1a2, 0, col);
			Meio1(pmeio1, fA1a2, A1a2, 1, col);
			Meio1(pmeio1, fA1a2, A1a2, 2, col);
			
		}
	
		///////////////
		//Meiosis II//
		/////////////
		CopyMat(aa,faa);
		//CopyMat(AA,fAA);
		//CopyMat(Aa,fAa);
		CopyMat(A1a,fA1a);
		//CopyMat(A1A,fA1A);
		//CopyMat(a2a,fa2a);
		//CopyMat(Aa2,fAa2);
		CopyMat(A1a2,fA1a2);
		
		for (col=0; col<(M/2+1); col++) {
			Meio2(pmeio2, faa, aagam, 0, col);
			Meio2(pmeio2, faa, aagam, 1, col);
			
		//	Meio2(pmeio2, fAA, AAgam, 0, col);
		//	Meio2(pmeio2, fAA, AAgam, 1, col);
		
		//	Meio2(pmeio2, fAa, Aagam, 0, col);
		//	Meio2(pmeio2, fAa, Aagam, 1, col);
			
			Meio2(pmeio2, fA1a, A1agam, 0, col);
			Meio2(pmeio2, fA1a, A1agam, 1, col);
			
		//	Meio2(pmeio2, fA1A, A1Agam, 0, col);
		//	Meio2(pmeio2, fA1A, A1Agam, 1, col);
			
		//	Meio2(pmeio2, fa2a, a2agam, 0, col);
		//	Meio2(pmeio2, fa2a, a2agam, 1, col);
			
		//	Meio2(pmeio2, fAa2, Aa2gam, 0, col);
		//	Meio2(pmeio2, fAa2, Aa2gam, 1, col);
			
			Meio2(pmeio2, fA1a2, A1a2gam, 0, col);
			Meio2(pmeio2, fA1a2, A1a2gam, 1, col);
		}
		
		//Define each gamete matrix
		GetGametesLarge(aagam, AAgam, Aagam, A1agam, A1Agam, a2agam, Aa2gam, A1a2gam, agamete,  Agamete, a2gamete,  A1gamete); 
		
		//Syngamy
		for (col=0; col<(M+1); col++) {
			syngamy2(0, col, agamete, agamete, aa);
			syngamy2(1, col, agamete, agamete, aa);
			syngamy2(2, col, agamete, agamete, aa);
			
		//	syngamy2(0, col, Agamete, Agamete, AA);
		//	syngamy2(1, col, Agamete, Agamete, AA);
		//	syngamy2(2, col, Agamete, Agamete, AA);
			
		//	syngamy1(0, col, Agamete, agamete, Aa);
		//	syngamy1(1, col, Agamete, agamete, Aa);
		//	syngamy1(2, col, Agamete, agamete, Aa);
			
			syngamy1(0, col, A1gamete, agamete, A1a);
			syngamy1(1, col, A1gamete, agamete, A1a);
			syngamy1(2, col, A1gamete, agamete, A1a);
			
		//	syngamy1(0, col, A1gamete, Agamete, A1A);
		//	syngamy1(1, col, A1gamete, Agamete, A1A);
		//	syngamy1(2, col, A1gamete, Agamete, A1A);
			
		//	syngamy2x2(0, col, a2gamete, agamete, a2a);
		//	syngamy2x2(1, col, a2gamete, agamete, a2a);
		//	syngamy2x2(2, col, a2gamete, agamete, a2a);
			
		//	syngamy1(0, col, Agamete, a2gamete, Aa2);
		///	syngamy1(1, col, Agamete, a2gamete, Aa2);
		//  syngamy1(2, col, Agamete, a2gamete, Aa2);
			
			syngamy1(0, col, A1gamete, a2gamete, A1a2);
			syngamy1(1, col, A1gamete, a2gamete, A1a2);
			syngamy1(2, col, A1gamete, a2gamete, A1a2);
		}

		//Normalize
		Normalize(aa, Aa, AA, A1a, A1A, a2a, Aa2, A1a2);			
	}
	
	//write data into file
	WriteMemoryFile(memory, 1, index);
	WriteMemoryFile(fitnessmem, 2, index);
//	WriteMemoryFile(variancemem, 3, index);
	
	//WriteGenotypeFile(aa, 1, index);
	//WriteGenotypeFile(Aa, 2, index);
	//WriteGenotypeFile(AA, 3, index);
	//WriteGenotypeFile(A1a, 4, index);
	//WriteGenotypeFile(A1A, 5, index);
	//WriteGenotypeFile(a2a, 6, index);
	//WriteGenotypeFile(Aa2, 7, index);
	//WriteGenotypeFile(A1a2, 8, index);
	
	//WriteHomFile(hom, index);

	return 0;	
}



