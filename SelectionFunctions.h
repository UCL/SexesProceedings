
 double  w(int nuc,  double  mit);

 double  wbar( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1]);

 double  wbarType( double  mat[3][M+1]);

void Select( double mata[3][M+1], double matb[3][M+1], double matc[3][M+1]);

double  wbarLarge( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1],  double  matd[3][M+1],  double  mate[3][M+1],  double  matf[3][M+1],  double  matg[3][M+1],  double  math[3][M+1]);

void SelectionLarge( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1],  double  matFixedd[3][M+1],  double  matChanged[3][M+1], double  matFixede[3][M+1],  double  matChangee[3][M+1], double  matFixedf[3][M+1],  double  matChangef[3][M+1], double  matFixedg[3][M+1],  double  matChangeg[3][M+1], double  matFixedh[3][M+1],  double  matChangeh[3][M+1]);

 double  wbarLargeRecombination( double  mat1[3][M+1],  double  mat2[3][M+1],  double  mat3[3][M+1],  double  mat4[3][M+1],  double  mat5[3][M+1],  double  mat6[3][M+1],  double  mat7[3][M+1],  double  mat8[3][M+1],  double  mat9[3][M+1],  double  mat10[3][M+1],  double  mat11[3][M+1],  double  mat12[3][M+1],  double  mat13[3][M+1],  double  mat14[3][M+1],  double  mat15[3][M+1]);

void SelectionRecombination( double  matFixed1[3][M+1],  double  matChange1[3][M+1],  double  matFixed2[3][M+1],  double  matChange2[3][M+1],  double  matFixed3[3][M+1],  double  matChange3[3][M+1],  double  matFixed4[3][M+1],  double  matChange4[3][M+1], double  matFixed5[3][M+1],  double  matChange5[3][M+1], double  matFixed6[3][M+1],  double  matChang6[3][M+1], double  matFixed7[3][M+1],  double  matChange7[3][M+1], double  matFixed8[3][M+1],  double  matChange8[3][M+1],  double  matFixed9[3][M+1],  double  matChange9[3][M+1],  double  matFixed10[3][M+1],  double  matChange10[3][M+1],  double  matFixed11[3][M+1],  double  matChange11[3][M+1],  double  matFixed12[3][M+1],  double  matChange12[3][M+1],  double  matFixed13[3][M+1],  double  matChange13[3][M+1],  double  matFixed14[3][M+1],  double  matChange14[3][M+1],  double  matFixed15[3][M+1],  double  matChange15[3][M+1]);

void SelectWithinType( double  matFixeda[3][M+1],  double  matChangea[3][M+1],  double  matFixedb[3][M+1],  double  matChangeb[3][M+1],  double  matFixedc[3][M+1],  double  matChangec[3][M+1]);

 double  wbarWithinType( double  mat[3][M+1]);

void CorrectMatrix(double matChange[3][M+1], double sum);