#include <R.h>
#include <Rdefines.h>
#include "catdap.h"

extern void F77_NAME(catdap1f)(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, int*, double*, int*);

SEXP catdap1m(SEXP nsamp, SEXP n, SEXP l, SEXP recode, SEXP n1, SEXP iskip1, SEXP n4, SEXP item1, SEXP item2, SEXP iconv, SEXP ires, SEXP iex, SEXP iskip, SEXP cdata)
{
    double *d1,*d2;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17,*i18;

    SEXP ans = R_NilValue, nc = R_NilValue, ia = R_NilValue, p = R_NilValue, total = R_NilValue, aic = R_NilValue, ord = R_NilValue;

    double *xp, *xaic = NULL;
    int  *xnc, *xia, *xtotal, *xord = NULL;

    int ll, nn, nn4, n4ln, ln1;
    int i;

    i1 = INTEGER_POINTER(nsamp);
    i2 = INTEGER_POINTER(n);
    i3 = INTEGER_POINTER(l);
    i4 = INTEGER_POINTER(recode);
    i5 = INTEGER_POINTER(n1);
    i6 = INTEGER_POINTER(iskip1);
    i7 = INTEGER_POINTER(n4);
    i8 = INTEGER_POINTER(item1);
    i9 = INTEGER_POINTER(item2);
    i10 = INTEGER_POINTER(iconv);
    i11 = INTEGER_POINTER(ires);
    i12 = INTEGER_POINTER(iex);
    i13 = INTEGER_POINTER(iskip);
    i14 = INTEGER_POINTER(cdata);

    nn = *i2;
    ll = *i3;
    nn4 = *i7;
    n4ln = nn4 * nn4 * ll * nn;
    ln1 = ll * (nn-1);

    PROTECT(ans = allocVector(VECSXP, 6));
    SET_VECTOR_ELT(ans, 0, nc = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 1, ia = allocVector(INTSXP, n4ln));
    SET_VECTOR_ELT(ans, 2, p = allocVector(REALSXP, n4ln));
    SET_VECTOR_ELT(ans, 3, total = allocVector(INTSXP, nn*nn4));
    SET_VECTOR_ELT(ans, 4, aic = allocVector(REALSXP, ll*nn));
    SET_VECTOR_ELT(ans, 5, ord = allocVector(INTSXP, ln1));

    i15 = INTEGER_POINTER(nc);
    i16 = INTEGER_POINTER(ia);
    d1 = NUMERIC_POINTER(p);
    i17 = INTEGER_POINTER(total);
    d2 = NUMERIC_POINTER(aic);
    i18 = INTEGER_POINTER(ord);

    F77_CALL(catdap1f) (i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,d1,i17,d2,i18);

    xnc = INTEGER(nc);
    xia = INTEGER(ia);
    xp = REAL(p);
    xtotal = INTEGER(total);
    xaic = REAL(aic);
    xord = INTEGER(ord);

    for(i=0; i<nn; i++) xnc[i] = i15[i];
    for(i=0; i<n4ln; i++) xia[i] = i16[i];
    for(i=0; i<n4ln; i++) xp[i] = d1[i];
    for(i=0; i<nn*nn4; i++) xtotal[i] = i17[i];
    for(i=0; i<ll*nn; i++) xaic[i] = d2[i];
    for(i=0; i<ln1; i++) xord[i] = i18[i];

    UNPROTECT(1);

    return ans;
}

