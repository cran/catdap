#include <R.h>
#include <Rdefines.h>
#include "catdap.h"

extern void F77_NAME(catdap2mf)(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, double*, int*, double*, int*, double*, double*, int*, int*, int*, int*, double*, int*, double*, int*, int*, double*, double*, int*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*, int*, int*, int*, int*, int*, double*, int*);

SEXP catdap2m(SEXP nsamp, SEXP n, SEXP l, SEXP recode, SEXP iskip1, SEXP it, SEXP nov, SEXP icl, SEXP item1, SEXP item2, SEXP pool, SEXP iconv, SEXP ires, SEXP iskip, SEXP isk, SEXP sk, SEXP icls, SEXP xx, SEXP ida, SEXP da, SEXP typeu, SEXP nmiss, SEXP n11, SEXP n33, SEXP ikr, SEXP jkr, SEXP ikkk, SEXP eps)
{
    double *d1,*d2,*d3,*d4,*d5,*d6,*d7,*d8,*d9,*d10,*d11,*d12;
    int *i1,*i2,*i3,*i4,*i5,*i6,*i7,*i8,*i9,*i10,*i11,*i12,*i13,*i14,*i15,*i16,*i17,*i18,*i19,*i20,*i21,*i22,*i23,*i24,*i25,*i26,*i27,*i28,*i29,*i30,*i31,*i32,*i33,*i34,*i35,*i36,*i37,*i38;

    SEXP ans = R_NilValue, iab = R_NilValue, totalc = R_NilValue, ttrr = R_NilValue, ab = R_NilValue, iaa = R_NilValue, pa = R_NilValue, idata = R_NilValue, ite = R_NilValue, dx = R_NilValue, aaam = R_NilValue, caa = R_NilValue, icaa = R_NilValue, nnia = R_NilValue, lk77 = R_NilValue, morder = R_NilValue, iby = R_NilValue, ibc = R_NilValue, pbc = R_NilValue, aic1 = R_NilValue, iabse = R_NilValue, baic = R_NilValue, ier = R_NilValue;

    double *xab, *xpa, *xdx, *xaaam, *xpbc, *xaic1, *xbaic = NULL;
    int  *xiab, *xtotalc, *xttrr, *xiaa, *xidata, *xite, *xcaa, *xicaa, *xnnia, *xlk77, *xmorder, *xiby, *xibc, *xiabse, *xier = NULL;

    int nn, nn1, nn3, nn31, nvm, iikr, jjkr, iikkk, icl1;
    int i;

    i1 = INTEGER_POINTER(nsamp);
    i2 = INTEGER_POINTER(n);
    i3 = INTEGER_POINTER(l);
    i4 = INTEGER_POINTER(recode);
    i5 = INTEGER_POINTER(iskip1);
    i6 = INTEGER_POINTER(it);
    i7 = INTEGER_POINTER(nov);
    i8 = INTEGER_POINTER(icl);
    i9 = INTEGER_POINTER(item1);
    i10 = INTEGER_POINTER(item2);
    i11 = INTEGER_POINTER(pool);
    i12 = INTEGER_POINTER(iconv);
    i13 = INTEGER_POINTER(ires);
    i14 = INTEGER_POINTER(iskip);
    i15 = INTEGER_POINTER(isk);
    d1 = NUMERIC_POINTER(sk);
    i16 = INTEGER_POINTER(icls);
    d2 = NUMERIC_POINTER(xx);
    i17 = INTEGER_POINTER(ida);
    d3 = NUMERIC_POINTER(da);
    d4 = NUMERIC_POINTER(typeu);
    i18 = INTEGER_POINTER(nmiss);
    i33 = INTEGER_POINTER(n11);
    i34 = INTEGER_POINTER(n33);
    i35 = INTEGER_POINTER(ikr);
    i36 = INTEGER_POINTER(jkr);
    i37 = INTEGER_POINTER(ikkk); 
    d12 = NUMERIC_POINTER(eps);

    nn = *i2;
    nvm = *i7;
    icl1 = *i8 + 1;
    nn1 = *i33;
    nn3 = *i34;
    iikr = *i35;
    jjkr = *i36;
    iikkk = *i37;
    nn31 = iikkk * icl1;

    PROTECT(ans = allocVector(VECSXP, 22));
    SET_VECTOR_ELT(ans, 0, iab = allocVector(INTSXP, nn*nn1*nn3));
    SET_VECTOR_ELT(ans, 1, totalc = allocVector(INTSXP, nn1));
    SET_VECTOR_ELT(ans, 2, ttrr = allocVector(INTSXP, nn*nn3));
    SET_VECTOR_ELT(ans, 3, ab = allocVector(REALSXP, nn*nn3));
    SET_VECTOR_ELT(ans, 4, iaa = allocVector(INTSXP, nn*nn3));
    SET_VECTOR_ELT(ans, 5, pa = allocVector(REALSXP, nn*nn1*nn3));
    SET_VECTOR_ELT(ans, 6, idata = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 7, ite = allocVector(INTSXP, nn));
    SET_VECTOR_ELT(ans, 8, dx = allocVector(REALSXP, nn));
    SET_VECTOR_ELT(ans, 9, aaam = allocVector(REALSXP, iikr));
    SET_VECTOR_ELT(ans, 10, caa = allocVector(INTSXP, iikr*jjkr));
    SET_VECTOR_ELT(ans, 11, icaa = allocVector(INTSXP, iikr));
    SET_VECTOR_ELT(ans, 12, nnia = allocVector(INTSXP, iikr));
    SET_VECTOR_ELT(ans, 13, lk77 = allocVector(INTSXP, 1));
    SET_VECTOR_ELT(ans, 14, morder = allocVector(INTSXP, iikr));
    SET_VECTOR_ELT(ans, 15, iby = allocVector(INTSXP, nvm*nn31));
    SET_VECTOR_ELT(ans, 16, ibc = allocVector(INTSXP, nn1*nn31));
    SET_VECTOR_ELT(ans, 17, pbc = allocVector(REALSXP, nn1*nn31));
    SET_VECTOR_ELT(ans, 18, aic1 = allocVector(REALSXP, icl1));
    SET_VECTOR_ELT(ans, 19, iabse = allocVector(INTSXP, nvm*nn3*icl1));
    SET_VECTOR_ELT(ans, 20, baic = allocVector(REALSXP, 1));
    SET_VECTOR_ELT(ans, 21, ier = allocVector(INTSXP, 3));

    i19 = INTEGER_POINTER(iab);
    i20 = INTEGER_POINTER(totalc);
    i21 = INTEGER_POINTER(ttrr);
    d5 = NUMERIC_POINTER(ab);
    i22 = INTEGER_POINTER(iaa);
    d6 = NUMERIC_POINTER(pa);
    i23 = INTEGER_POINTER(idata);
    i24 = INTEGER_POINTER(ite);
    d7 = NUMERIC_POINTER(dx);
    d8 = NUMERIC_POINTER(aaam);
    i25 = INTEGER_POINTER(caa);
    i26 = INTEGER_POINTER(icaa);
    i27 = INTEGER_POINTER(nnia);
    i28 = INTEGER_POINTER(lk77);
    i29 = INTEGER_POINTER(morder);
    i30 = INTEGER_POINTER(iby);
    i31 = INTEGER_POINTER(ibc);
    d9 = NUMERIC_POINTER(pbc);
    d10 = NUMERIC_POINTER(aic1);
    i32 = INTEGER_POINTER(iabse);
    d11 = NUMERIC_POINTER(baic);
    i38 = INTEGER_POINTER(ier);

    F77_CALL(catdap2mf) (i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,d1,i16,d2,
i17,d3,d4,i18,i19,i20,i21,d5,i22,d6,i23,i24,d7,d8,i25,i26,i27,i28,i29,i30,i31,d9,d10,i32,d11,i33,i34,i35,i36,i37,d12,i38);

    xiab = INTEGER(iab);
    xtotalc = INTEGER(totalc);
    xttrr = INTEGER(ttrr);
    xab = REAL(ab);
    xiaa = INTEGER(iaa);
    xpa = REAL(pa);
    xidata = INTEGER(idata);
    xite = INTEGER(ite);
    xdx = REAL(dx);
    xaaam = REAL(aaam);
    xcaa = INTEGER(caa);
    xicaa = INTEGER(icaa);
    xnnia = INTEGER(nnia);
    xlk77 = INTEGER(lk77);
    xmorder = INTEGER(morder);
    xiby = INTEGER(iby);
    xibc = INTEGER(ibc);
    xpbc = REAL(pbc);
    xaic1 = REAL(aic1);
    xiabse = INTEGER(iabse);
    xbaic = REAL(baic);
    xier = INTEGER(ier);

    for(i=0; i<nn*nn1*nn3; i++) xiab[i] = i19[i];
    for(i=0; i<nn1; i++) xtotalc[i] = i20[i];
    for(i=0; i<nn*nn3; i++) xttrr[i] = i21[i];
    for(i=0; i<nn*nn3; i++) xab[i] = d5[i];
    for(i=0; i<nn*nn3; i++) xiaa[i] = i22[i];
    for(i=0; i<nn*nn1*nn3; i++) xpa[i] = d6[i];
    for(i=0; i<nn; i++) xidata[i] = i23[i];
    for(i=0; i<nn; i++) xite[i] = i24[i];
    for(i=0; i<nn; i++) xdx[i] = d7[i];
    for(i=0; i<iikr; i++) xaaam[i] = d8[i];
    for(i=0; i<iikr*jjkr; i++) xcaa[i] = i25[i];
    for(i=0; i<iikr; i++) xicaa[i] = i26[i];
    for(i=0; i<iikr; i++) xnnia[i] = i27[i];
    *xlk77 = *i28;
    for(i=0; i<iikr; i++) xmorder[i] = i29[i];
    for(i=0; i<nvm*nn31; i++) xiby[i] = i30[i];
    for(i=0; i<nn1*nn31; i++) xibc[i] = i31[i];
    for(i=0; i<nn1*nn31; i++) xpbc[i] = d9[i];
    for(i=0; i<icl1; i++) xaic1[i] = d10[i];
    for(i=0; i<nvm*nn3*icl1; i++) xiabse[i] = i32[i];
    *xbaic = *d11;
    for(i=0; i<3; i++) xier[i] = i38[i];

    UNPROTECT(1);

    return ans;
}
