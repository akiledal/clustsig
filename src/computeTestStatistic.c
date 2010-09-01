/* computeTestStatistic.c */

#include "computeTestStatistic.h"

SEXP computeTestStatistic(SEXP data, SEXP averages, SEXP length) {
  double *Pdata, *Paverages, *teststat; 
  int p = 0, len; // PROTECT tracker
  SEXP ts; 

  // this should be a 2xlength array of sorts
  PROTECT(data = AS_NUMERIC(data)); p++;
  Pdata = NUMERIC_POINTER(data);

  PROTECT(averages = AS_NUMERIC(averages)); p++;
  Paverages = NUMERIC_POINTER(averages);

  PROTECT(ts = NEW_NUMERIC(1)); p++;
  teststat = NUMERIC_POINTER(ts);
  teststat[0] = 0;

  len = INTEGER_VALUE(length);

  for (int i = 0; i < len; i++){
    void R_CheckUserInterrupt(void);
    teststat[0] += fabs(Pdata[0+2*i] - Paverages[i]);
  }

  UNPROTECT(p);
  return ts;
}
