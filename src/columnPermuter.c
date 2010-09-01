#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

SEXP columnPermuter(SEXP, SEXP, SEXP);
void shuffle(int, double *, double *, int);

SEXP columnPermuter(SEXP matrix, SEXP nrow, SEXP ncol){
  double *mat, *newmat;
  int n = 0;
  int row = INTEGER_VALUE(nrow);
  int col = INTEGER_VALUE(ncol);
  SEXP newmatrix;

  PROTECT(matrix = AS_NUMERIC(matrix)); n++;
  mat = NUMERIC_POINTER(matrix);

  PROTECT(newmatrix = allocMatrix(REALSXP, row, col)); n++;
  newmat = REAL(newmatrix);
  
  double temp[row];
  for (int j = 0; j < col; j++){
    for (int i = 0; i < row; i++){
      void R_CheckUserInterrupt(void);
      temp[i]=mat[i+row*j];
    }
    shuffle(row, newmat, temp, j); // newmat receives the shuffle (passed by reference) 
  }

  UNPROTECT(n);
  return(newmatrix);
}

void shuffle(int length, double *res, double *temp, int curcol){
  int x = length;
  int j;
  for (int i = 0; i < length; i++) {
    GetRNGstate();
    j = x * unif_rand();
    PutRNGstate();
    res[i + curcol*length] = temp[j];
    temp[j] = temp[--x];
  }
}
