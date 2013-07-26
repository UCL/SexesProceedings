#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "SelectionFunctions.h"

////////////////////////////
//FUNTIONS FOR GENERAL USE//
///////////////////////////

//function for a choose b
 double  com( double  a,  double  b)
{
	 double  result;
	int i;
	 double  p;
	if(a<b || b<0)
		result=0;
	else if(b==0)
		result=1;
	else{
		p=1.0;
		for(i=0; i<b; i++)
			p=p*((a-i)/(b-i));
		result = p;
	}
	return 	result;
}//com

//row sumation function
 double  SumRow( double  *A, int ncols, int row){
	 double  sum=0.0;
	int COL;
	
	for (A +=ncols*row, COL = 0; COL<ncols; COL++) sum+=*A++;
	return(sum);
}//SumRow

//column sumation function
 double  SumCol( double  *A, int nrows, int ncols, int col){
	 double  sum=0.0;
	int ROW;
	
	for (A += col, ROW=0; ROW<nrows; A += ncols, ROW++) sum+= *A; 
	
	return(sum);
}// SumCol

//maximum of two integers
int max(int a, int b){
	if (a>b)return(a);
	else return(b);
}//max

//minimum of two integers
int min(int a, int b){
	if(a<b)return(a);
	else return(b);
}//min

//copy matrix a onto matrix b
void CopyMat( double  a[3][M+1],  double  b[3][M+1]){
	int row;
	int col;
	for (row=0; row<3; row++) 
		for (col=0; col<(M+1); col++)
			b[row][col] = a[row][col];
}//CopyMat

//copy matrix a onto matrix b for larger matrices - use in development step
void CopyMatLarge( double  a[3][2*M+1],  double  b[3][2*M+1]){
	int row;
	int col;
	for (row=0; row<3; row++) 
		for (col=0; col<(2*M+1); col++)
			b[row][col] = a[row][col];
}//CopyMatLarge

//function to sum over a matrix
double  SumMatrix(double  mat[3][M+1]){
	int row;
	int col;
	double  sum=0.0;
	
	for (row=0; row<3; row++) 
		for (col=0; col<(M+1); col++) sum += mat[row][col];
}//SumMatrix

//function to sum over a smaller matrix
 double  SumMatrixSmall( double  mat[2][M+1]){
	int row;
	int col;
	 double  sum=0.0;
	
	for (row=0; row<2; row++) 
		for (col=0; col<(M+1); col++) sum += mat[row][col];
}//SumMatrixSmall

//Normalize matrices after syngmay for small problem
void NormalizeSmall( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1]){
	double  sum = 0.0;
	int row, col;
	
	sum = SumMatrix(mata) + SumMatrix(matb) + SumMatrix(matc);
	
	for(row=0; row<3; row++){
		for(col=0;col<(M+1);col++){
			mata[row][col]= mata[row][col]/sum;
			matb[row][col]= matb[row][col]/sum;
			matc[row][col]= matc[row][col]/sum;
		}
	}
}//Normalize



//Normalize matrices at the end
void Normalize( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1],  double  matd[3][M+1],  double  mate[3][M+1],  double  matf[3][M+1],  double  matg[3][M+1],  double  math[3][M+1]){
	 double  sum = 0.0;
	int row, col;
	
	sum = SumMatrix(mata) + SumMatrix(matb) + SumMatrix(matc) + SumMatrix(matd) + SumMatrix(mate) + SumMatrix(matf) + SumMatrix(matg) + SumMatrix(math);
	
	for(row=0; row<3; row++){
		for(col=0;col<(M+1);col++){
			mata[row][col]= mata[row][col]/sum;
			matb[row][col]= matb[row][col]/sum;
			matc[row][col]= matc[row][col]/sum;
			matd[row][col]= matd[row][col]/sum;
			mate[row][col]= mate[row][col]/sum;
			matf[row][col]= matf[row][col]/sum;
			matg[row][col]= matg[row][col]/sum;
			math[row][col]= math[row][col]/sum;
		}
	}
}//Normalize

//Normalize matrices at the end in the recombination model
void NormalizeRecombination( double  mat1[3][M+1],  double  mat2[3][M+1],  double  mat3[3][M+1],  double  mat4[3][M+1],  double  mat5[3][M+1],  double  mat6[3][M+1],  double  mat7[3][M+1],  double  mat8[3][M+1],  double  mat9[3][M+1],  double  mat10[3][M+1],  double  mat11[3][M+1],  double  mat12[3][M+1],  double  mat13[3][M+1],  double  mat14[3][M+1],  double  mat15[3][M+1]){
	 double  sum = 0.0;
	int row, col;
	
	sum = SumMatrix(mat1) + SumMatrix(mat2) + SumMatrix(mat3) + SumMatrix(mat4) + SumMatrix(mat5) + SumMatrix(mat6) + SumMatrix(mat7) + SumMatrix(mat8)+ SumMatrix(mat9) + SumMatrix(mat10) + SumMatrix(mat11) + SumMatrix(mat12) + SumMatrix(mat13) + SumMatrix(mat14) + SumMatrix(mat15);
	
	for(row=0; row<3; row++){
		for(col=0;col<(M+1);col++){
			mat1[row][col]= mat1[row][col]/sum;
			mat2[row][col]= mat2[row][col]/sum;
			mat3[row][col]= mat3[row][col]/sum;
			mat4[row][col]= mat4[row][col]/sum;
			mat5[row][col]= mat5[row][col]/sum;
			mat6[row][col]= mat6[row][col]/sum;
			mat7[row][col]= mat7[row][col]/sum;
			mat8[row][col]= mat8[row][col]/sum;
			mat9[row][col]= mat9[row][col]/sum;
			mat10[row][col]= mat10[row][col]/sum;
			mat11[row][col]= mat11[row][col]/sum;
			mat12[row][col]= mat12[row][col]/sum;
			mat13[row][col]= mat13[row][col]/sum;
			mat14[row][col]= mat14[row][col]/sum;
			mat15[row][col]= mat15[row][col]/sum;
		}
	}
}//NormalizeRecombination

//Finding the variance of a type 
 double  VarType( double  mat[3][M+1]){
	 double  meansq=0.0;
	 double  variance = 0.0;
	 double  mean = 0.0;
	int col;
	
	for (col=0; col<(M+1); col++) meansq += mat[0][col]*pow(w(0, col),2) + mat[1][col]*pow(w(1,col),2) + mat[2][col]*pow(w(2,col),2);
	meansq = meansq/(SumMatrix(mat));
	mean = wbarType(mat);
	variance = meansq - pow(mean,2);
	
	return(variance);
}//VarType
	
//sum over a col of the basic matrices
 double  SumNuc( double  mat[3][M+1], int col){
	 double  sum=0.0;
	int row;
	
	for (row=0; row<(M+1); row++) sum += mat[col][row];
	
	return(sum);
}
//SumNuc
	
//maximum frequency of Aa in memory
int MemMax( double  mat[runs][genotypes], int row){
	 double  max = 0;
	int i, maxindex=0;
	
	for (i=0; i<runs; i++){ 
		if(mat[i][row]>max){
			max=mat[i][row];
			maxindex = i;
		}
	}	
	return(maxindex);
}//MemMaxAa

//minimum frequency of Aa in memory
int MemMin( double  mat[runs][genotypes], int row){
	 double  min = 10;
	int i, minindex=0;
	
	for (i=1500; i<(runs-10); i++){ 
		if(mat[i][row]<min) 
		{
			min=mat[i][row];
			minindex = i;
		}
	}
	return(minindex);
}//MemMaxAa



