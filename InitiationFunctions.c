#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include "Global.h"
#include "GeneralFunctions.h"
#include "InitiationFunctions.h"


///////////////////////////////////
/////////Initiating matirces//////
/////////////////////////////////
//General matrix initialization
void InitMatrix(double  Table[NUC][MIT+1]){
	int nuc=0;
	int mit=0;
	for(nuc=0; nuc<3; nuc++)
		for(mit=0;mit<(MIT+1);mit++)
			Table[nuc][mit]=0;
}//InitMatrix

//Initialize matrix a
void InitMatrixa(double  Table[NUC][MIT+1]){
	int nuc=0;
	int mit=0;
	double  sum = 0.0;
	
	for(nuc=0; nuc<1; nuc++)
		for(mit=0;mit<(M+1);mit++){
			//Table[nuc][mit]=1.0; 
			Table[nuc][mit]=rand(); 
			//std::cout << "rand: " << rand() << std::endl;
			sum += Table[nuc][mit];
		}
	for(nuc=0; nuc<1; nuc++)
		for(mit=0;mit<(M+1);mit++){
			Table[nuc][mit]=Table[nuc][mit]/sum;
		}
	
}//InitMatrixa

//Initialize matrix b
void InitMatrixb( double  Table[NUC][MIT+1]){
	int nuc=0;
	int mit=0;
	for(nuc=0; nuc<NUC; nuc++)
		for(mit=0;mit<(MIT+1);mit++)
			Table[nuc][mit]=1.0/(24.0*(M+1));
}//InitMatrixb

//Initialize matrix c
void InitMatrixc( double  Tablea[NUC][MIT+1],  double  Tablec[NUC][MIT+1],  double  ratio){
	int nuc=0;
	int mit=0;
	for(nuc=0; nuc<NUC; nuc++)
		for(mit=0;mit<(MIT+1);mit++){
			Tablec[nuc][mit]= Tablea[nuc][mit]/ratio;
			//Tablea[nuc][mit]= Tablea[nuc][mit] - Tablea[nuc][mit]/ratio;
			Tablea[nuc][mit]= Tablea[nuc][mit] - Tablec[nuc][mit];
		}
}//InitMatrixc

//Initialize gamete matrices
void InitGameteMatrix( double  Table[2][MIT+1]){
	int nuc=0;
	int mit=0;
	for(nuc=0; nuc<2; nuc++)
		for(mit=0;mit<(MIT+1);mit++)
			Table[nuc][mit]=0;
}//InitMatrix

//Initialize memory matrix
void InitMem( double  Table[runs][3]){
	int row=0;
	int col=0;
	for(row=0; row<(runs); row++)
		for(col=0;col<3;col++)
			Table[row][col]=0;
}//InitMem

//Initialize large memory matrix
void InitMemLarge( double  Table[runs][8]){
	int row=0;
	int col=0;
	for(row=0; row<(runs); row++)
		for(col=0;col<8;col++)
			Table[row][col]=0;
}//InitMemLarge

//Initialize homozygote state matrix
void InitHom( double  Table[runs][1]){
	int row=0;
	for(row=0; row<(runs); row++)
			Table[row][0]=0;
}//InitHom


//Initialize array for developmental states
void InitArray( double  array[3*(M+1)][3][(M+1)]){
	int i,j,k;
	
	for (i=0; i<3*(M+1); i++) 
		for (j=0; j<3; j++) 
			for (k=0; k<(M+1); k++) 
				array[i][j][k]=0;
}
	
//Initialize larger matrix
void InitMatrixLarge( double  Table[3][2*M+1]){
	int nuc=0;
	int mit=0;
	for(nuc=0; nuc<3; nuc++)
		for(mit=0;mit<(2*M+1);mit++)
			Table[nuc][mit]=0;
}//InitMatrix
	
	
	
	

