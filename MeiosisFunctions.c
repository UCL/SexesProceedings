#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "MutationFunctions.h"
#include "SelectionFunctions.h"
#include "MeiosisFunctions.h"


///////////////////////////////////
///Funstions relevant to MEIOSIS//
/////////////////////////////////
//initiating the meiosis step 1 matrix 
void InitMeio1( double  mat[M+1][M+1]){
	int row, col;
	for(row=0; row<(M+1); row++)
		for (col=0; col<(M+1); col++) 
			mat[row][col]=com(2.0*row,col)*com(2.0*(M-row), M-col)/com(2.0*M,M);
}//InitMeio1

//initiating the meiosis step 2 matrix 
void InitMeio2( double  mat[M+1][M+1]){
	int row, col;
	for(row=0; row<(M+1); row++)
		for (col=0; col<(M+1); col++) 
			mat[row][col]=com(row, col)*com(M-row, M/2.0-col)/com(M,M/2.0);
}//InitMeio2

//Tramformation following the first meiotic step
void Meio1( double  meioprobs1[M+1][M+1],  double  matFixed[3][M+1],  double  matChange[3][M+1], int nuc, int mit){
	 double  sum1 = 0.0;
	 double  sum2 = 0.0;
	 double  sum = 0.0;
	int k;
	
	switch (nuc) {
		case (0):
			for (k=0; k<(M+1); k++) {
				sum1 += meioprobs1[k][mit]*matFixed[nuc][k];
				sum2 += meioprobs1[k][mit]*matFixed[1][k]*(1.0/6.0);
			}
			matChange[nuc][mit] = sum1+sum2;
			break;
		case (2):
			for (k=0; k<(M+1); k++) {
				sum1 += meioprobs1[k][mit]*matFixed[nuc][k];
				sum2 += meioprobs1[k][mit]*matFixed[1][k]*(1.0/6.0);
			}
			matChange[nuc][mit] = sum1+sum2;
			break;
		default:
			for (k=0; k<(M+1); k++) sum+= meioprobs1[k][mit]*matFixed[nuc][k]*(2.0/3.0);
			matChange[nuc][mit] = sum;
			break;
	}
}//Meio1

//Tramformation following the second meiotic step
void Meio2( double  meioprobs2[M+1][M+1],  double  matFixed[3][M+1],  double  matChange[2][M+1] ,int nuc, int mit){
	 double  sum1 = 0.0;
	 double  sum2 = 0.0;
	int k; //index for mitochondria summation 
	
	switch (nuc) {
		case (0):
			for (k=0; k<(M+1); k++) {
				sum1 += meioprobs2[k][mit]*matFixed[0][k];
				sum2 += meioprobs2[k][mit]*matFixed[1][k]*0.5;
			}
			matChange[nuc][mit] = sum1+sum2;
			break;
		default:
			sum1=0.0;
			sum2=0.0;
			
			for (k=0; k<(M+1); k++){
				sum1 += meioprobs2[k][mit]*matFixed[2][k];
				sum2 += meioprobs2[k][mit]*matFixed[1][k]*0.5;
			}
			matChange[nuc][mit] = sum1+sum2;
			break;
	}
}//Meio2



