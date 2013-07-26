#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "MutationFunctions.h"
#include "BiogenesisFunctions.h"

//Functions here are used to impose selfish mutants. This is essentially a mitochondrial
//biogenesis step where mutant mitochondria are given a multiplicative advantage. 


//probability function for biogenesis
 double  pbio( double  x){
	 double  k=0.5;
	return((1.0+k)*x/(M+x*k));
}

//initialize biogenesis probability matrix
void InitBioMat( double  Table[M+1][M+1]){
	int row, col;
	
	for (row=0; row<(M+1); row++) {
		for (col=0; col<(M+1); col++) {
			Table[row][col] = com(M, col)*pow(pbio(row), col)*pow(1.0 - pbio(row), M-col);
		}
	}
	for (col=0; col<M; col++) {
		Table[M][col] = 0.0;
	}
	Table[M][M]=1.0;
	Table[0][0]=1.0;
}

//biogenesis function
void Biogenesis( double  matfixed[3][M+1], double  matchange[3][M+1], double  biomat[M+1][M+1]){
	 double  sum=0.0;
	int j, nuc, mit;
	
	for (nuc=0; nuc<1; nuc++) {
		for (mit=0; mit<(M+1); mit++){
			sum=0.0;
			for (j=0; j<(M+1); j++) {
				sum += biomat[j][mit]*matfixed[nuc][j];
			}
			matchange[nuc][mit] = sum;
		}
	}		
}
