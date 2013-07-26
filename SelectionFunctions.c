#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "MutationFunctions.h"
#include "SelectionFunctions.h"

////////////////////////////////////
//Funstions relevant to SELECTION//
//////////////////////////////////
//fitness function
 double  w(int nuc,  double  mit){
	
	 double  tolerance=1.0;
	 double  p = (mit)/( double )M;
	p = p/tolerance;
	switch (nuc) {
		case 0:
			return(1 - pow(p,0.5));		
			break;
		case 1:
			//if(p>0.5) return(1 - pow(p, 2.0));
			//else return(1 - pow(1-p, 2.0));
			return(1 - 0.5*pow(p,0.5) - 0.5*pow(1-p,0.5));
		default:
			return(1 - pow(1-p,0.5));		
			break;
	}//switch
}//w


//get mean fitness
 double  wbar( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1]){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mata[0][col]*w(0, col) + mata[1][col]*w(1,col) + mata[2][col]*w(2,col)+matb[0][col]*w(0, col) + matb[1][col]*w(1,col) + matb[2][col]*w(2,col)+matc[0][col]*w(0, col) + matc[1][col]*w(1,col) + matc[2][col]*w(2,col);
	return(sum);
}//wbar

//get mean fitness for larger problem
 double  wbarLarge( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1],  double  matd[3][M+1],  double  mate[3][M+1],  double  matf[3][M+1],  double  matg[3][M+1],  double  math[3][M+1]){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mata[0][col]*w(0, col) + mata[1][col]*w(1,col) + mata[2][col]*w(2,col)+matb[0][col]*w(0, col) + matb[1][col]*w(1,col) + matb[2][col]*w(2,col)+matc[0][col]*w(0, col) + matc[1][col]*w(1,col) + matc[2][col]*w(2,col)+matd[0][col]*w(0, col) + matd[1][col]*w(1,col) + matd[2][col]*w(2,col)+mate[0][col]*w(0, col) + mate[1][col]*w(1,col) + mate[2][col]*w(2,col) +matf[0][col]*w(0, col) + matf[1][col]*w(1,col) + matf[2][col]*w(2,col)+matg[0][col]*w(0, col) + matg[1][col]*w(1,col) + matg[2][col]*w(2,col)+math[0][col]*w(0, col) + math[1][col]*w(1,col) + math[2][col]*w(2,col);
	return(sum);
}//wbarLarge

//get mean fitness for larger problem with recombination
 double  wbarLargeRecombination( double  mat1[3][M+1],  double  mat2[3][M+1],  double  mat3[3][M+1],  double  mat4[3][M+1],  double  mat5[3][M+1],  double  mat6[3][M+1],  double  mat7[3][M+1],  double  mat8[3][M+1],  double  mat9[3][M+1],  double  mat10[3][M+1],  double  mat11[3][M+1],  double  mat12[3][M+1],  double  mat13[3][M+1],  double  mat14[3][M+1],  double  mat15[3][M+1]){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mat1[0][col]*w(0, col) + mat1[1][col]*w(1,col) + mat1[2][col]*w(2,col)+mat2[0][col]*w(0, col) + mat2[1][col]*w(1,col) + mat2[2][col]*w(2,col)+mat3[0][col]*w(0, col) + mat3[1][col]*w(1,col) + mat3[2][col]*w(2,col)+mat4[0][col]*w(0, col) + mat4[1][col]*w(1,col) + mat4[2][col]*w(2,col)+mat5[0][col]*w(0, col) + mat5[1][col]*w(1,col) + mat5[2][col]*w(2,col) +mat6[0][col]*w(0, col) + mat6[1][col]*w(1,col) + mat6[2][col]*w(2,col)+mat7[0][col]*w(0, col) + mat7[1][col]*w(1,col) + mat7[2][col]*w(2,col)+mat8[0][col]*w(0, col) + mat8[1][col]*w(1,col) + mat8[2][col]*w(2,col) + mat9[0][col]*w(0, col) + mat9[1][col]*w(1,col) + mat9[2][col]*w(2,col)+mat10[0][col]*w(0, col) + mat10[1][col]*w(1,col) + mat10[2][col]*w(2,col)+mat11[0][col]*w(0, col) + mat11[1][col]*w(1,col) + mat11[2][col]*w(2,col)+mat12[0][col]*w(0, col) + mat12[1][col]*w(1,col) + mat12[2][col]*w(2,col) +mat13[0][col]*w(0, col) + mat13[1][col]*w(1,col) + mat13[2][col]*w(2,col)+mat14[0][col]*w(0, col) + mat14[1][col]*w(1,col) + mat14[2][col]*w(2,col)+mat15[0][col]*w(0, col) + mat15[1][col]*w(1,col) + mat15[2][col]*w(2,col);
	return(sum);
}//wbarLargeRecombination

//get mean fitness for one type
 double  wbarType( double  mat[3][M+1]){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mat[0][col]*w(0, col) + mat[1][col]*w(1,col) + mat[2][col]*w(2,col);
	return(sum/SumMatrix(mat));		
}//wbarType

//apply selection
void Select( double mata[3][M+1], double matb[3][M+1], double matc[3][M+1]){
	 double  mean = 1.0;
	int col;
	
	mean = wbar(mata, matb, matc);
	for (col=0; col<(M+1); col++) {
		mata[0][col] = mata[0][col]*w(0, col)/mean;
		mata[1][col] = mata[1][col]*w(1, col)/mean;
		mata[2][col] = mata[2][col]*w(2, col)/mean;
		
		matb[0][col] = matb[0][col]*w(0, col)/mean;
		matb[1][col] = matb[1][col]*w(1, col)/mean;
		matb[2][col] = matb[2][col]*w(2, col)/mean;
		
		matc[0][col] = matc[0][col]*w(0, col)/mean;
		matc[1][col] = matc[1][col]*w(1, col)/mean;
		matc[2][col] = matc[2][col]*w(2, col)/mean;
		
	}//for
}//Selection

//apply selection to larger problem
void SelectionLarge( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1],  double  matFixedd[3][M+1],  double  matChanged[3][M+1], double  matFixede[3][M+1],  double  matChangee[3][M+1], double  matFixedf[3][M+1],  double  matChangef[3][M+1], double  matFixedg[3][M+1],  double  matChangeg[3][M+1], double  matFixedh[3][M+1],  double  matChangeh[3][M+1]){
	 double  mean = 1.0;
	int col;
	
	mean = wbarLarge(matFixeda, matFixedb, matFixedc, matFixedd, matFixede, matFixedf, matFixedg, matFixedh);
	for (col=0; col<(M+1); col++) {
		matChangea[0][col] = matFixeda[0][col]*w(0, col)/mean;
		matChangea[1][col] = matFixeda[1][col]*w(1, col)/mean;
		matChangea[2][col] = matFixeda[2][col]*w(2, col)/mean;
		
		matChangeb[0][col] = matFixedb[0][col]*w(0, col)/mean;
		matChangeb[1][col] = matFixedb[1][col]*w(1, col)/mean;
		matChangeb[2][col] = matFixedb[2][col]*w(2, col)/mean;
		
		matChangec[0][col] = matFixedc[0][col]*w(0, col)/mean;
		matChangec[1][col] = matFixedc[1][col]*w(1, col)/mean;
		matChangec[2][col] = matFixedc[2][col]*w(2, col)/mean;
		
		matChanged[0][col] = matFixedd[0][col]*w(0, col)/mean;
		matChanged[1][col] = matFixedd[1][col]*w(1, col)/mean;
		matChanged[2][col] = matFixedd[2][col]*w(2, col)/mean;
		
		matChangee[0][col] = matFixede[0][col]*w(0, col)/mean;
		matChangee[1][col] = matFixede[1][col]*w(1, col)/mean;
		matChangee[2][col] = matFixede[2][col]*w(2, col)/mean;
		matChangef[0][col] = matFixedf[0][col]*w(0, col)/mean;
		matChangef[1][col] = matFixedf[1][col]*w(1, col)/mean;
		matChangef[2][col] = matFixedf[2][col]*w(2, col)/mean;
		
		matChangeg[0][col] = matFixedg[0][col]*w(0, col)/mean;
		matChangeg[1][col] = matFixedg[1][col]*w(1, col)/mean;
		matChangeg[2][col] = matFixedg[2][col]*w(2, col)/mean;
		
		matChangeh[0][col] = matFixedh[0][col]*w(0, col)/mean;
		matChangeh[1][col] = matFixedh[1][col]*w(1, col)/mean;
		matChangeh[2][col] = matFixedh[2][col]*w(2, col)/mean;
		
	}//for
}//SelectionLarge

//apply selection to model with recombination
void SelectionRecombination( double  matFixed1[3][M+1],  double  matChange1[3][M+1],  double  matFixed2[3][M+1],  double  matChange2[3][M+1],  double  matFixed3[3][M+1],  double  matChange3[3][M+1],  double  matFixed4[3][M+1],  double  matChange4[3][M+1], double  matFixed5[3][M+1],  double  matChange5[3][M+1], double  matFixed6[3][M+1],  double  matChange6[3][M+1], double  matFixed7[3][M+1],  double  matChange7[3][M+1], double  matFixed8[3][M+1],  double  matChange8[3][M+1],  double  matFixed9[3][M+1],  double  matChange9[3][M+1],  double  matFixed10[3][M+1],  double  matChange10[3][M+1],  double  matFixed11[3][M+1],  double  matChange11[3][M+1],  double  matFixed12[3][M+1],  double  matChange12[3][M+1],  double  matFixed13[3][M+1],  double  matChange13[3][M+1],  double  matFixed14[3][M+1],  double  matChange14[3][M+1],  double  matFixed15[3][M+1],  double  matChange15[3][M+1]){
	 double  mean = 1.0;
	int col;
	
	mean = wbarLargeRecombination(matFixed1, matFixed2, matFixed3, matFixed4, matFixed5, matFixed6, matFixed7, matFixed8, matFixed9, matFixed10, matFixed11, matFixed12, matFixed13, matFixed14, matFixed15);
	for (col=0; col<(M+1); col++) {
		matChange1[0][col] = matFixed1[0][col]*w(0, col)/mean;
		matChange1[1][col] = matFixed1[1][col]*w(1, col)/mean;
		matChange1[2][col] = matFixed1[2][col]*w(2, col)/mean;
		
		matChange2[0][col] = matFixed2[0][col]*w(0, col)/mean;
		matChange2[1][col] = matFixed2[1][col]*w(1, col)/mean;
		matChange2[2][col] = matFixed2[2][col]*w(2, col)/mean;
		
		matChange3[0][col] = matFixed3[0][col]*w(0, col)/mean;
		matChange3[1][col] = matFixed3[1][col]*w(1, col)/mean;
		matChange3[2][col] = matFixed3[2][col]*w(2, col)/mean;
		
		matChange4[0][col] = matFixed4[0][col]*w(0, col)/mean;
		matChange4[1][col] = matFixed4[1][col]*w(1, col)/mean;
		matChange4[2][col] = matFixed4[2][col]*w(2, col)/mean;
		
		matChange5[0][col] = matFixed5[0][col]*w(0, col)/mean;
		matChange5[1][col] = matFixed5[1][col]*w(1, col)/mean;
		matChange5[2][col] = matFixed5[2][col]*w(2, col)/mean;
		
		matChange6[0][col] = matFixed6[0][col]*w(0, col)/mean;
		matChange6[1][col] = matFixed6[1][col]*w(1, col)/mean;
		matChange6[2][col] = matFixed6[2][col]*w(2, col)/mean;
		
		matChange7[0][col] = matFixed7[0][col]*w(0, col)/mean;
		matChange7[1][col] = matFixed7[1][col]*w(1, col)/mean;
		matChange7[2][col] = matFixed7[2][col]*w(2, col)/mean;
		
		matChange8[0][col] = matFixed8[0][col]*w(0, col)/mean;
		matChange8[1][col] = matFixed8[1][col]*w(1, col)/mean;
		matChange8[2][col] = matFixed8[2][col]*w(2, col)/mean;
		
		matChange9[0][col] = matFixed9[0][col]*w(0, col)/mean;
		matChange9[1][col] = matFixed9[1][col]*w(1, col)/mean;
		matChange9[2][col] = matFixed9[2][col]*w(2, col)/mean;
		
		matChange10[0][col] = matFixed10[0][col]*w(0, col)/mean;
		matChange10[1][col] = matFixed10[1][col]*w(1, col)/mean;
		matChange10[2][col] = matFixed10[2][col]*w(2, col)/mean;
		
		matChange11[0][col] = matFixed11[0][col]*w(0, col)/mean;
		matChange11[1][col] = matFixed11[1][col]*w(1, col)/mean;
		matChange11[2][col] = matFixed11[2][col]*w(2, col)/mean;
		
		matChange12[0][col] = matFixed12[0][col]*w(0, col)/mean;
		matChange12[1][col] = matFixed12[1][col]*w(1, col)/mean;
		matChange12[2][col] = matFixed12[2][col]*w(2, col)/mean;
		
		matChange13[0][col] = matFixed13[0][col]*w(0, col)/mean;
		matChange13[1][col] = matFixed13[1][col]*w(1, col)/mean;
		matChange13[2][col] = matFixed13[2][col]*w(2, col)/mean;
		
		matChange14[0][col] = matFixed14[0][col]*w(0, col)/mean;
		matChange14[1][col] = matFixed14[1][col]*w(1, col)/mean;
		matChange14[2][col] = matFixed14[2][col]*w(2, col)/mean;
		
		matChange15[0][col] = matFixed15[0][col]*w(0, col)/mean;
		matChange15[1][col] = matFixed15[1][col]*w(1, col)/mean;
		matChange15[2][col] = matFixed15[2][col]*w(2, col)/mean;
		
	}//for
}//SelectionRecombination

//FOR THE MULTICELLULAR CASE
//weight function
 double  weight(int nuc, int mit){
	
	switch (nuc) {
		case 0:
			return(exp(-mit/(M+0.0)));
			break;
		case 2:
			return(exp(-(M-mit)/(M+0.0)));
			break;
		default:
			return(0.5*exp(-mit/(M+0.0)) + 0.5*exp(-(M-mit)/(M+0.0)));
			break;
	}
}//weight

//get multicellular mean
 double  MultMean( double  mat[3][M+1]){
	 double  sum=0.0;
	 double  mean=0.0;
	int i,j;
	
	for (i=0; i<3; i++) 
		for (j=0; j<(M+1); j++) sum += mat[i][j]*weight(i,j);	

	for (i=0; i<3; i++) 
		for (j=0; j<(M+1); j++) mean += w(i,j)*mat[i][j]*weight(i, j)/sum;
	
	return(mean);
}//NormWeights


//the functions below will be used to impose selection withou changing the frequency of each genotype in the simple case 
//for the purposes of getting t plot usefull for the paper
//selection over each type

//get mean fitness for one type
 double  wbarWithinType( double  mat[3][M+1]){
	int col;
	 double  sum=0.0;
	
	for (col=0; col<(M+1); col++) sum += mat[0][col]*w(0, col) + mat[1][col]*w(1,col) + mat[2][col]*w(2,col);
	return(sum);		
}//wbarType

//apply selection
void SelectWithinType( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1]){
	 double  mean1 = 1.0;
	 double  mean2 = 1.0;
	 double  mean3 = 1.0;
	int col;
	
	mean1 = wbarType(matFixeda);
	mean2 = wbarType(matFixedb);
	mean3 = wbarType(matFixedc);
	
	//mean1 = wbarWithinType(matFixeda);
	//mean2 = wbarWithinType(matFixedb);
	//mean3 = wbarWithinType(matFixedc);
	
	for (col=0; col<(M+1); col++) {
		matChangea[0][col] = matFixeda[0][col]*w(0, col)/mean1;
		matChangea[1][col] = matFixeda[1][col]*w(1, col)/mean1;
		matChangea[2][col] = matFixeda[2][col]*w(2, col)/mean1;
		
		matChangeb[0][col] = matFixedb[0][col]*w(0, col)/mean2;
		matChangeb[1][col] = matFixedb[1][col]*w(1, col)/mean2;
		matChangeb[2][col] = matFixedb[2][col]*w(2, col)/mean2;
		
		matChangec[0][col] = matFixedc[0][col]*w(0, col)/mean3;
		matChangec[1][col] = matFixedc[1][col]*w(1, col)/mean3;
		matChangec[2][col] = matFixedc[2][col]*w(2, col)/mean3;
	}//for
}//Selection

void CorrectMatrix(double matChange[3][M+1], double sum){
	double total=0.0;
	double error = 0.0;
	int col;
	
	total = SumMatrix(matChange);
	error = (sum/total);
	for (col=0; col<(M+1); col++) {
		matChange[0][col] = matChange[0][col]*error;
	}
}

		
	
	




