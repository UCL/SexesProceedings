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
#include "SyngamyFunctions.h"



///////////////////////////////////
///Funstions relevant to SYNGAMY//
/////////////////////////////////
//extracting gamete matrices
void GetGametes( double  mata[2][M+1],  double  matb[2][M+1],  double  matc[2][M+1],  double  areturn[2][M+1],  double  Areturn[2][M+1]){
	int row, col;
	
	for (row=0; row<2; row++) {
		for (col=0; col<(M+1); col++) {
			areturn[row][col] = mata[row][col] + 0.5*matc[row][col];
			Areturn[row][col] = matb[row][col] + 0.5*matc[row][col];
 		}
	}//for
}//GetGametes

//extracting gamete matrices for larger problem
void GetGametesLarge( double  mata[2][M+1],  double  matb[2][M+1],  double  matc[2][M+1],  double  matd[2][M+1],  double  mate[2][M+1],  double  matf[2][M+1], double  matg[2][M+1],  double  math[2][M+1],  double  areturn[2][M+1],  double  Areturn[2][M+1], double  a2return[2][M+1],  double  A1return[2][M+1]){
	int row, col;
	
	for (row=0; row<2; row++) {
		for (col=0; col<(M+1); col++) {
			areturn[row][col] = mata[row][col] + 0.5*matc[row][col] + 0.5*matd[row][col] + 0.5*matf[row][col];
			Areturn[row][col] = matb[row][col] + 0.5*matc[row][col] + 0.5*mate[row][col] + 0.5*matg[row][col];
			a2return[row][col] = 0.5*matf[row][col] + 0.5*matg[row][col] + 0.5*math[row][col];
			A1return[row][col] = 0.5*matd[row][col] + 0.5*mate[row][col] + 0.5*math[row][col];
 		}
	}//for
}//GetGametesLarge

//extracting gamete matrices for larger problem with recombination
void GetGametesRecombination( double  matA1A2[2][M+1],  double  matA1a2[2][M+1],  double  matA1a[2][M+1],  double  matA1A[2][M+1],  double  matA2a1[2][M+1],  double  mata1a2[2][M+1],  double  mata1a[2][M+1],  double  matAa1[2][M+1],  double  matA2a[2][M+1],  double  matA2A[2][M+1],  double  mata2a[2][M+1],  double  matAa2[2][M+1],  double  matAa[2][M+1],  double  mataa[2][M+1],  double  matAA[2][M+1],  double  areturn[2][M+1],  double  Areturn[2][M+1], double  a2return[2][M+1],  double  A2return[2][M+1],  double  a1return[2][M+1],  double  A1return[2][M+1]){
	int row, col;
	
	for (row=0; row<2; row++) {
		for (col=0; col<(M+1); col++) {
			areturn[row][col] = (1.0/4.0)*matA1a[row][col] + (1.0/2.0)*mata1a[row][col] + (1.0/4.0)*matAa1[row][col] + (1.0/4.0)*matA2a[row][col] + (1.0/2.0)*mata2a[row][col] + (1.0/4.0)*matAa2[row][col] + (1.0/2.0)*matAa[row][col] + mataa[row][col];
			Areturn[row][col] = (1.0/4.0)*matA1a[row][col] + (1.0/2.0)*matA1A[row][col] + (1.0/4.0)*matAa1[row][col] + (1.0/4.0)*matA2a[row][col] + (1.0/2.0)*matA2A[row][col] + (1.0/4.0)*matAa2[row][col] + (1.0/2.0)*matAa[row][col] + matAA[row][col];
			a1return[row][col] = (1.0/4.0)*matA1a2[row][col] + (1.0/4.0)*matA1a[row][col] + (1.0/4.0)*matA2a1[row][col] + (1.0/2.0)*mata1a2[row][col] + (1.0/2.0)*mata1a[row][col] + (1.0/4.0)*matAa1[row][col];
			A1return[row][col] = (1.0/2.0)*matA1A2[row][col] + (1.0/4.0)*matA1a2[row][col] + (1.0/4.0)*matA1a[row][col] + (1.0/2.0)*matA1A[row][col] + (1.0/4.0)*matA2a1[row][col] + (1.0/4.0)*matAa1[row][col];
			a2return[row][col] = (1.0/4.0)*matA1a2[row][col] + (1.0/4.0)*matA2a1[row][col] + (1.0/2.0)*mata1a2[row][col] + (1.0/4.0)*matA2a[row][col] + (1.0/2.0)*mata2a[row][col] + (1.0/4.0)*matAa2[row][col];
			A2return[row][col] = (1.0/2.0)*matA1A2[row][col] + (1.0/4.0)*matA1a2[row][col] + (1.0/4.0)*matA2a1[row][col] + (1.0/4.0)*matA2a[row][col] + (1.0/2.0)*matA2A[row][col] + (1.0/4.0)*matAa2[row][col];
 		}
	}//for
}//GetGametesRecombination 

//syngamy function - uniparental
void syngamy1(int nuc,int mit, double  act[2][M+1], double  pass[2][M+1],  double  returnmat[3][M+1]){
	 double  sum=0.0;
	 double  sum1=0.0;
	 double  sum2=0.0;
	 double  coef=0.0;
	 double  coef1=0.0;
	 double  coef2=0.0;
	int l;
	
	switch (nuc) {
		case 0:
			coef = SumRow(&pass[0][0], M+1, 0);
			for (l=0; l<(M/2+1); l++) sum += coef*act[0][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			returnmat[nuc][mit] = 2.0*sum;
			break;
			
		case 1:
			coef1 = SumRow(&pass[0][0], M+1, 0);
			coef2 = SumRow(&pass[0][0], M+1, 1);
			for (l=0; l<(M/2+1); l++){ 
				sum1 += coef1*act[1][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
				sum2 += coef2*act[0][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			}
			returnmat[nuc][mit] = 2.0*(sum1+sum2);
			break;
		default:
			coef = SumRow(&pass[0][0], M+1, 1);
			for (l=0; l<(M/2+1); l++) sum += coef*act[1][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			returnmat[nuc][mit] = 2.0*sum;
			break;
	}//switch
}//syngamy1

//syngamy function - biparental
void syngamy2(int nuc,int mit, double  act1[2][M+1], double  act2[2][M+1],  double  returnmat[3][M+1]){
	 double  sum=0.0;
	int l;
	
	switch (nuc) {
		case 0:
			for (l=0; l<(mit+1); l++) sum += act1[0][l]*act2[0][mit-l];
			returnmat[nuc][mit] = sum;
			break;
			
		case 1:
			for (l=0; l<(mit+1); l++) sum += act1[0][l]*act2[1][mit-l] + act2[0][l]*act1[1][mit-l];
			returnmat[nuc][mit] = sum;
			break;
		default:
			for (l=0; l<(mit+1); l++) sum += act1[1][l]*act2[1][mit-l];
			returnmat[nuc][mit] =  sum;
			break;
	}//switch
}//syngamy2

//syngamy function - biparental with a multiple of 2 because of different gametes 
void syngamy2x2(int nuc,int mit, double  act1[2][M+1], double  act2[2][M+1],  double  returnmat[3][M+1]){
	 double  sum=0.0;
	int l;
	
	switch (nuc) {
		case 0:
			for (l=0; l<(mit+1); l++) sum += act1[0][l]*act2[0][mit-l];
			returnmat[nuc][mit] = 2*sum;
			break;
			
		case 1:
			for (l=0; l<(mit+1); l++) sum += act1[0][l]*act2[1][mit-l] + act2[0][l]*act1[1][mit-l];
			returnmat[nuc][mit] = 2*sum;
			break;
		default:
			for (l=0; l<(mit+1); l++) sum += act1[1][l]*act2[1][mit-l];
			returnmat[nuc][mit] =  2*sum;
			break;
	}//switch
}//syngamy2x2

//syngamy function - uniparental but with two A gametes so no need to multiply by 2
void syngamy1AxA(int nuc,int mit, double  act[2][M+1], double  pass[2][M+1],  double  returnmat[3][M+1]){
	double  sum=0.0;
	double  sum1=0.0;
	double  sum2=0.0;
	double  coef=0.0;
	double  coef1=0.0;
	double  coef2=0.0;
	int l;
	
	switch (nuc) {
		case 0:
			coef = SumRow(&pass[0][0], M+1, 0);
			for (l=0; l<(M/2+1); l++) sum += coef*act[0][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			returnmat[nuc][mit] = sum;
			break;
			
		case 1:
			coef1 = SumRow(&pass[0][0], M+1, 0);
			coef2 = SumRow(&pass[0][0], M+1, 1);
			for (l=0; l<(M/2+1); l++){ 
				sum1 += coef1*act[1][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
				sum2 += coef2*act[0][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			}
			returnmat[nuc][mit] = (sum1+sum2);
			break;
		default:
			coef = SumRow(&pass[0][0], M+1, 1);
			for (l=0; l<(M/2+1); l++) sum += coef*act[1][l]*pow((l/(M/2.0)), mit)*pow(((M/2.0) - l)/(M/2.0), (M-mit))*com(M, mit);
			returnmat[nuc][mit] = sum;
			break;
	}//switch
}//syngamy1


