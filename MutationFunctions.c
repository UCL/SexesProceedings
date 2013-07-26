#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "MutationFunctions.h"


///////////////////////////////////
//Funstions relevant to MUTATION//
/////////////////////////////////
//Binomial probability density function
 double  BinomDens( double  x,  double  N,  double  prob){
	return(com(N,x)*pow(prob, x)*pow(1-prob,N-x));
}//BinomDens
//Mitochondrial Mutation Summation
 double  SumMitProbs(int l, int j){
	 double  s=0.0;
	int x;
	for (x=max(0, j-l); x<(min(M-l,j)+1); x++) s += BinomDens(x, M-l,mu)*BinomDens(l-j+x, l, mub);
	return(s);
}//SumMitProbs

//initiating the mitochondrial mutation matrix 
void InitMutMat( double  mat[M+1][M+1]){
	int l,j;
	for (l=0; l<(M+1); l++) 
		for (j=0; j<(M+1); j++) 
			mat[l][j]=SumMitProbs(l,j);
}//InitMutMat

//Nuclear mutation rate from i to j
 double  NucMut(int i, int j){
	switch (i) {
		case 0:
			switch (j) {
				case 0:
					return(pow(1-nu, 2));
					break;
				case 1:
					return(2.0*(1-nu)*nu);
					break;
				default:
					return(pow(nu,2));
					break;
			}
			break;
		case 1:
			switch (j) {
				case 0:
					return((1-nu)*nub);
					break;
				case 1:
					return((1-nu)*(1-nub)+nub*nu);
					break;
				default:
					return((1-nub)*nu);
					break;
			}
			break;
		default:
			switch (j) {
				case 0:
					return(pow(nub,2));
					break;
				case 1:
					return(2.0*(1-nub)*nub);
					break;
				default:
					return(pow(1-nub,2));
					break;
			}
			break;
	}//switch   
}//NucMut

//Mutation function - overall
void Mutation( double  probs[M+1][M+1], double  matFixed[3][M+1], double  matChange[3][M+1], int i, int j){
	 double  mut = 0.0;
	int k;
	
	for (k=0; k<(M+1); k++) mut += NucMut(0,i)*probs[k][j]*matFixed[0][k] + NucMut(1,i)*probs[k][j]*matFixed[1][k] + NucMut(2,i)*probs[k][j]*matFixed[2][k];
	
	matChange[i][j] = mut;
	
}//Mutation

