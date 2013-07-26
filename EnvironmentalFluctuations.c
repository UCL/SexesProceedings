#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"

///////////////////////////////////////////////////////////////////////////////////////////
//////////////FUNTIONS FOR IMPLEMENTING SELECTION UNDER A CHANGING ENVIRONMENT////////////
/////////////////////////////////////////////////////////////////////////////////////////


//These functions are used to impose a changing environment in the mitonuclear co-adaptation framwork.
//References to small and large problem correspond to the simple aa+A case versus aa+a1+A2 case respectively. 

//fitness function with environmental fluctuations
 double  wenv(int nuc,  double  mit, int index,  double  homozygote){
	
	 double  tolerance=1.0;
	 double  result = 0.0;
	 double  p = (mit)/( double )M;
	p = p/tolerance;
	switch (nuc) {
		case 0:
			result = 1 - pow(p,2.0);
			if (index%changenuc<duration && index>300) 
				if (homozygote>0.5) {
					if (result-E>0.0){
						result = result-E;
					}					
					else {
						result = 0.0;
					}
				}
			return(result);
			break;
			
		case 1:
			return(1 - 0.5*pow(p,2.0) - 0.5*pow(1-p,2.0));
			break;
			
		default:
			result = 1 - pow(1-p,2.0);
			if (index%changenuc<duration && index>300) 
				if (homozygote<0.5) {
					if (result-E>0.0){
						result = result-E;
					}
					else {
						result = 0.0;
					}
				}
			return(result);
			break;
	}//switch
}//wenv


//get mean fitness for smaller problem considering environmental changes
 double  wbarEnv( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1], int index,  double  homozygote){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mata[0][col]*wenv(0, col, index, homozygote) + mata[1][col]*wenv(1,col, index, homozygote) + mata[2][col]*wenv(2,col, index, homozygote)+matb[0][col]*wenv(0, col, index, homozygote) + matb[1][col]*wenv(1,col, index, homozygote) + matb[2][col]*wenv(2,col, index, homozygote)+matc[0][col]*wenv(0, col, index, homozygote) + matc[1][col]*wenv(1,col, index, homozygote) + matc[2][col]*wenv(2,col, index, homozygote);
	return(sum);
}//wbarEnv

//apply selection to smaller problem considering environmental changes
void SelectEnv( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1], int index,  double  homozygote){
	 double  mean = 1.0;
	int col;
	
	mean = wbarEnv(matFixeda, matFixedb, matFixedc, index, homozygote);
	for (col=0; col<(M+1); col++) {
		matChangea[0][col] = matFixeda[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangea[1][col] = matFixeda[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangea[2][col] = matFixeda[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangeb[0][col] = matFixedb[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangeb[1][col] = matFixedb[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangeb[2][col] = matFixedb[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangec[0][col] = matFixedc[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangec[1][col] = matFixedc[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangec[2][col] = matFixedc[2][col]*wenv(2, col, index, homozygote)/mean;
	}//for
}//SelectEnv


//get mean fitness for larger problem considering environmental changes
 double  wbarLargeEnv( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1],  double  matd[3][M+1],  double  mate[3][M+1],  double  matf[3][M+1],  double  matg[3][M+1],  double  math[3][M+1], int index,  double  homozygote){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mata[0][col]*wenv(0, col, index, homozygote) + mata[1][col]*wenv(1,col, index, homozygote) + mata[2][col]*wenv(2,col, index, homozygote)+matb[0][col]*wenv(0, col, index, homozygote) + matb[1][col]*wenv(1,col, index, homozygote) + matb[2][col]*wenv(2,col, index, homozygote)+matc[0][col]*wenv(0, col, index, homozygote) + matc[1][col]*wenv(1,col, index, homozygote) + matc[2][col]*wenv(2,col, index, homozygote)+matd[0][col]*wenv(0, col, index, homozygote) + matd[1][col]*wenv(1,col, index, homozygote) + matd[2][col]*wenv(2,col, index, homozygote)+mate[0][col]*wenv(0, col, index, homozygote) + mate[1][col]*wenv(1,col, index, homozygote) + mate[2][col]*wenv(2,col, index, homozygote)+matf[0][col]*wenv(0, col, index, homozygote) + matf[1][col]*wenv(1,col, index, homozygote) + matf[2][col]*wenv(2,col, index, homozygote)+matg[0][col]*wenv(0, col, index, homozygote) + matg[1][col]*wenv(1,col, index, homozygote) + matg[2][col]*wenv(2,col, index, homozygote)+math[0][col]*wenv(0, col, index, homozygote) + math[1][col]*wenv(1,col, index, homozygote) + math[2][col]*wenv(2,col, index, homozygote);
	return(sum);
}//wbarLargeEnv

//get mean fitness for larger problem with recombination and environmental fluctuations
 double  wbarRecombinationEnv( double  mat1[3][M+1],  double  mat2[3][M+1],  double  mat3[3][M+1],  double  mat4[3][M+1],  double  mat5[3][M+1],  double  mat6[3][M+1],  double  mat7[3][M+1],  double  mat8[3][M+1],  double  mat9[3][M+1],  double  mat10[3][M+1],  double  mat11[3][M+1],  double  mat12[3][M+1],  double  mat13[3][M+1],  double  mat14[3][M+1],  double  mat15[3][M+1], int index,  double  homozygote){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mat1[0][col]*wenv(0, col, index, homozygote) + mat1[1][col]*wenv(1,col, index, homozygote) + mat1[2][col]*wenv(2,col, index, homozygote)+mat2[0][col]*wenv(0, col, index, homozygote) + mat2[1][col]*wenv(1,col, index, homozygote) + mat2[2][col]*wenv(2,col, index, homozygote)+mat3[0][col]*wenv(0, col, index, homozygote) + mat3[1][col]*wenv(1,col, index, homozygote) + mat3[2][col]*wenv(2,col, index, homozygote)+mat4[0][col]*wenv(0, col, index, homozygote) + mat4[1][col]*wenv(1,col, index, homozygote) + mat4[2][col]*wenv(2,col, index, homozygote)+mat5[0][col]*wenv(0, col, index, homozygote) + mat5[1][col]*wenv(1,col, index, homozygote) + mat5[2][col]*wenv(2,col, index, homozygote) +mat6[0][col]*wenv(0, col, index, homozygote) + mat6[1][col]*wenv(1,col, index, homozygote) + mat6[2][col]*wenv(2,col, index, homozygote)+mat7[0][col]*wenv(0, col, index, homozygote) + mat7[1][col]*wenv(1,col, index, homozygote) + mat7[2][col]*wenv(2,col, index, homozygote)+mat8[0][col]*wenv(0, col, index, homozygote) + mat8[1][col]*wenv(1,col, index, homozygote) + mat8[2][col]*wenv(2,col, index, homozygote) + mat9[0][col]*wenv(0, col, index, homozygote) + mat9[1][col]*wenv(1,col, index, homozygote) + mat9[2][col]*wenv(2,col, index, homozygote)+mat10[0][col]*wenv(0, col, index, homozygote) + mat10[1][col]*wenv(1,col, index, homozygote) + mat10[2][col]*wenv(2,col, index, homozygote)+mat11[0][col]*wenv(0, col, index, homozygote) + mat11[1][col]*wenv(1,col, index, homozygote) + mat11[2][col]*wenv(2,col, index, homozygote)+mat12[0][col]*wenv(0, col, index, homozygote) + mat12[1][col]*wenv(1,col, index, homozygote) + mat12[2][col]*wenv(2,col, index, homozygote) +mat13[0][col]*wenv(0, col, index, homozygote) + mat13[1][col]*wenv(1,col, index, homozygote) + mat13[2][col]*wenv(2,col, index, homozygote)+mat14[0][col]*wenv(0, col, index, homozygote) + mat14[1][col]*wenv(1,col, index, homozygote) + mat14[2][col]*wenv(2,col, index, homozygote)+mat15[0][col]*wenv(0, col, index, homozygote) + mat15[1][col]*wenv(1,col, index, homozygote) + mat15[2][col]*wenv(2,col, index, homozygote);
	return(sum);
}//wbarRecombinationEnv

//apply selection to larger problem considering environmental changes
void SelectionLargeEnv( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1],  double  matFixedd[3][M+1],  double  matChanged[3][M+1], double  matFixede[3][M+1],  double  matChangee[3][M+1], double  matFixedf[3][M+1],  double  matChangef[3][M+1], double  matFixedg[3][M+1],  double  matChangeg[3][M+1], double  matFixedh[3][M+1],  double  matChangeh[3][M+1], int index,  double  homozygote){
	 double  mean = 1.0;
	int col;
	
	mean = wbarLargeEnv(matFixeda, matFixedb, matFixedc, matFixedd, matFixede, matFixedf, matFixedg, matFixedh, index, homozygote);
	for (col=0; col<(M+1); col++) {
		matChangea[0][col] = matFixeda[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangea[1][col] = matFixeda[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangea[2][col] = matFixeda[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangeb[0][col] = matFixedb[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangeb[1][col] = matFixedb[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangeb[2][col] = matFixedb[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangec[0][col] = matFixedc[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangec[1][col] = matFixedc[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangec[2][col] = matFixedc[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChanged[0][col] = matFixedd[0][col]*wenv(0, col, index, homozygote)/mean;
		matChanged[1][col] = matFixedd[1][col]*wenv(1, col, index, homozygote)/mean;
		matChanged[2][col] = matFixedd[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangee[0][col] = matFixede[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangee[1][col] = matFixede[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangee[2][col] = matFixede[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangef[0][col] = matFixedf[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangef[1][col] = matFixedf[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangef[2][col] = matFixedf[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangeg[0][col] = matFixedg[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangeg[1][col] = matFixedg[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangeg[2][col] = matFixedg[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChangeh[0][col] = matFixedh[0][col]*wenv(0, col, index, homozygote)/mean;
		matChangeh[1][col] = matFixedh[1][col]*wenv(1, col, index, homozygote)/mean;
		matChangeh[2][col] = matFixedh[2][col]*wenv(2, col, index, homozygote)/mean;
		
	}//for
}//SelectionLargeEnv

//apply selection to model with recombination and environmental fluctuations
void SelectionRecombinationEnv( double  matFixed1[3][M+1],  double  matChange1[3][M+1],  double  matFixed2[3][M+1],  double  matChange2[3][M+1],  double  matFixed3[3][M+1],  double  matChange3[3][M+1],  double  matFixed4[3][M+1],  double  matChange4[3][M+1], double  matFixed5[3][M+1],  double  matChange5[3][M+1], double  matFixed6[3][M+1],  double  matChange6[3][M+1], double  matFixed7[3][M+1],  double  matChange7[3][M+1], double  matFixed8[3][M+1],  double  matChange8[3][M+1],  double  matFixed9[3][M+1],  double  matChange9[3][M+1],  double  matFixed10[3][M+1],  double  matChange10[3][M+1],  double  matFixed11[3][M+1],  double  matChange11[3][M+1],  double  matFixed12[3][M+1],  double  matChange12[3][M+1],  double  matFixed13[3][M+1],  double  matChange13[3][M+1],  double  matFixed14[3][M+1],  double  matChange14[3][M+1],  double  matFixed15[3][M+1],  double  matChange15[3][M+1], int index,  double  homozygote){
	 double  mean = 1.0;
	int col;
	
	mean = wbarRecombinationEnv(matFixed1, matFixed2, matFixed3, matFixed4, matFixed5, matFixed6, matFixed7, matFixed8, matFixed9, matFixed10, matFixed11, matFixed12, matFixed13, matFixed14, matFixed15,index, homozygote);
	for (col=0; col<(M+1); col++) {
		matChange1[0][col] = matFixed1[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange1[1][col] = matFixed1[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange1[2][col] = matFixed1[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange2[0][col] = matFixed2[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange2[1][col] = matFixed2[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange2[2][col] = matFixed2[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange3[0][col] = matFixed3[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange3[1][col] = matFixed3[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange3[2][col] = matFixed3[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange4[0][col] = matFixed4[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange4[1][col] = matFixed4[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange4[2][col] = matFixed4[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange5[0][col] = matFixed5[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange5[1][col] = matFixed5[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange5[2][col] = matFixed5[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange6[0][col] = matFixed6[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange6[1][col] = matFixed6[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange6[2][col] = matFixed6[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange7[0][col] = matFixed7[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange7[1][col] = matFixed7[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange7[2][col] = matFixed7[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange8[0][col] = matFixed8[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange8[1][col] = matFixed8[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange8[2][col] = matFixed8[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange9[0][col] = matFixed9[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange9[1][col] = matFixed9[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange9[2][col] = matFixed9[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange10[0][col] = matFixed10[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange10[1][col] = matFixed10[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange10[2][col] = matFixed10[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange11[0][col] = matFixed11[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange11[1][col] = matFixed11[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange11[2][col] = matFixed11[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange12[0][col] = matFixed12[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange12[1][col] = matFixed12[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange12[2][col] = matFixed12[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange13[0][col] = matFixed13[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange13[1][col] = matFixed13[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange13[2][col] = matFixed13[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange14[0][col] = matFixed14[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange14[1][col] = matFixed14[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange14[2][col] = matFixed14[2][col]*wenv(2, col, index, homozygote)/mean;
		
		matChange15[0][col] = matFixed15[0][col]*wenv(0, col, index, homozygote)/mean;
		matChange15[1][col] = matFixed15[1][col]*wenv(1, col, index, homozygote)/mean;
		matChange15[2][col] = matFixed15[2][col]*wenv(2, col, index, homozygote)/mean;
		
	}//for
}//SelectionRecombinationEnv




//get mean fitness for one type during environmental fluctuations process
 double  wEnvbarType( double  mat[3][M+1], int index,  double  homozygote){
	int col;
	 double  sum=0.0;
	for (col=0; col<(M+1); col++) sum += mat[0][col]*wenv(0, col, index, homozygote) + mat[1][col]*wenv(1, col, index, homozygote) + mat[2][col]*wenv(2, col, index, homozygote);
	return(sum/SumMatrix(mat));	
}//wbarType




