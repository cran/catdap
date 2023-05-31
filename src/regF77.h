/*
*  Catdap : Categorical Data Analysis Program Package
*  Copyright (C) 2008    The Institute of Statistical Mathematics
*
*    This program is free software; you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation; either version 2 of the License, or
*    (at your option) any later version.
*
*    This program is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with this program; if not, write to the Free Software
*    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.
*
*  ismrp at grp.ism.ac.jp
*/

#include <R.h>
#include <Rinternals.h>
#include <libintl.h>

#define _(String) (String)

/* Fortran : */

void F77_NAME(catdap1)(int *nsamp, int *n, int *l, int *record, int *n1,
              int *iskip1, int *n4, int *item1, int *item2, int *iconv,
			  int *ires, int *iex, int *iskip, int *cdata, int *nc, int *ia,
              double *p, int *total, double *aic, int *ord);

void F77_NAME(catdap2m)(int *nsamp, int *n, int *l, int *record, int *iskip1,
              int *it, int *nov, int *icl, int *item1, int *item2, int *pool,
              int *iconv, int *ires, int *iskip, int *isk, double *sk,
              int *icls, double *xx, int *ida, double *da, double *typeu,
              int *nmiss, int *iab, int *totalc, int *ttrr, double *ab,
              int *iaa, double *pa, int *idata, int *ite, double *dx,
              double *aaam, int *caa, int *icaa, int *nnia, int *lk77,
              int *morder, int *ibc, double *pbc, double *aic1, int *iabse,
              double *baic, int *n11, int *n33, int *ikr, int *jkr, int *ikkk,
              double *eps, int *nrange, int *ier);
