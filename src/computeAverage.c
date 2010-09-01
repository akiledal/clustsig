#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>

SEXP computeAverage(SEXP expSP, SEXP numexp){
  double *Paverage, tempsum, *Pcurrent;
  SEXP average, current;
  int p = 0, num = INTEGER_VALUE(numexp), len=LENGTH(VECTOR_ELT(expSP, 1))/2;

  PROTECT(average = allocMatrix(REALSXP, 1, len)); p++;
  Paverage = REAL(average);
  
  for (int i = 0; i < len; i++){
    tempsum = 0;
    for (int j = 0; j < num; j++){
      void R_CheckUserInterrupt(void);
      PROTECT(current = AS_NUMERIC(VECTOR_ELT(expSP, j)));
      Pcurrent = REAL(current);
      tempsum += Pcurrent[2*i];
      UNPROTECT(1);
    }
    Paverage[i] = tempsum/num;
  }
  UNPROTECT(p);
  return(average);
}
