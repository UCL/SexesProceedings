
 double  BinomDens( double  x,  double  N,  double  prob);

 double  SumMitProbs(int l, int j);

void InitMutMat( double  mat[M+1][M+1]);

 double  NucMut(int i, int j);

void Mutation( double  probs[M+1][M+1], double  matFixed[3][M+1], double  matChange[3][M+1], int i, int j);