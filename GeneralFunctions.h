
 double  com( double  a,  double  b);

 double  SumRow( double  *A, int ncols, int row);

 double  SumCol( double  *A, int nrows, int ncols, int col);

int max(int a, int b);

int min(int a, int b);

void CopyMat( double  a[3][M+1],  double  b[3][M+1]);

 double  SumMatrix( double  mat[3][M+1]);

 double  SumMatrixSmall( double  mat[2][M+1]);

void NormalizeSmall( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1]);

void Normalize( double  mata[3][M+1],  double  matb[3][M+1],  double  matc[3][M+1],  double  matd[3][M+1],  double  mate[3][M+1],  double  matf[3][M+1],  double  matg[3][M+1],  double  math[3][M+1]);

void NormalizeRecombination( double  mat1[3][M+1],  double  mat2[3][M+1],  double  mat3[3][M+1],  double  mat4[3][M+1],  double  mat5[3][M+1],  double  mat6[3][M+1],  double  mat7[3][M+1],  double  mat8[3][M+1],  double  mat9[3][M+1],  double  mat10[3][M+1],  double  mat11[3][M+1],  double  mat12[3][M+1],  double  mat13[3][M+1],  double  mat14[3][M+1],  double  mat15[3][M+1]);

 double  VarType( double  mat[3][M+1]);

 double  SumNuc( double  mat[3][M+1], int col);

void CopyMatLarge( double  a[3][2*M+1],  double  b[3][2*M+1]);

int MemMax( double  mat[runs][genotypes], int row);

int MemMin( double  mat[runs][genotypes], int row);


