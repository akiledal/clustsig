/* braycurtis.c */
/* Progress: Mostly working. Need to implement data checking:
   Data should be type-checked
   Data should not have Missing/Na/NaN values. Figure out a behavior for undefined results (force a similarity of 1? give a logical option for users to choose?)
*/ 

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <math.h>

SEXP braycurtis(SEXP data, SEXP nrow, SEXP ncol, SEXP cval){
  double *Pdata, *Pdist;
  int p = 0; // to keep track of how many PROTECTs are used
  int row = INTEGER_VALUE(nrow);
  int col = INTEGER_VALUE(ncol);
  double c = NUMERIC_VALUE(cval);
  SEXP dist;
  
  PROTECT(data = AS_NUMERIC(data)); p++;
  Pdata = NUMERIC_POINTER(data);
  PROTECT(dist = allocMatrix(REALSXP, row, row)); p++; 
  Pdist = REAL(dist);

  int numerator, denominator;
  // i and j tell us which two rows/sites we are comparing
  for (int i = 0; i < row; i++){ 
    for (int j = i + 1; j < row; j++){
      numerator  = 0;
      denominator= 0;
      for (int k = 0; k < col; k++){
	void R_CheckUserInterrupt(void);
       	numerator += fabs(Pdata[i+k*row] - Pdata[j+k*row]);
	denominator += Pdata[i+k*row] + Pdata[j+k*row];
      }
      denominator += 2*c; // add in the dummy species, whatever it is
      if (denominator == 0) // in case c = 0
	Pdist[j+i*row] = 0;
      else
	Pdist[j+i*row] = 100.0 * numerator/denominator;
    }
  }
  UNPROTECT(p);
  return(dist);
}
