// ===============================================================
//
// SCArray R package
// Copyright (C) 2021-2023    Xiuwen Zheng
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "sc_main.h"


// ========================

static inline bool is_int(SEXP x)
{
	int v = TYPEOF(x);
	return (v==INTSXP) || (v==LGLSXP);
}

// SparseArraySeed
SparseMatrix::SparseMatrix(SEXP mat)
{
	static const char *err_type =
		"SparseArraySeed should be a numeric matrix.";
	// slot nzindex
	SEXP ii = GET_SLOT(mat, mkString("nzindex"));
	if (!Rf_isMatrix(ii)) Rf_error(err_type);
	nnzero = INTEGER(GET_DIM(ii))[0];
	nzi_r = INTEGER(ii); nzi_c = nzi_r + nnzero;
	// slot nzdata
	nzdata = GET_SLOT(mat, mkString("nzdata"));
	if (TYPEOF(nzdata)!=REALSXP && !is_int(nzdata))
		Rf_error(err_type);
	if (Rf_length(nzdata) != nnzero)
		Rf_error(err_type);
}


static void throw_error_type(SEXP x)
{
	SEXP nm = Rf_getAttrib(x, R_ClassSymbol);
	if (Rf_isString(nm))
	{
		const char *s = CHAR(STRING_ELT(nm, 0));
		Rf_error("Unsupported data type 'class: %s'.", s);
	} else {
		Rf_error("Unsupported data type 'typeof: %d'.", (int)TYPEOF(x));
	}
}


bool is_sparse_seed(SEXP mat)
{
	return Rf_inherits(mat, "SparseArraySeed") == TRUE;
}

void get_mat_size(SEXP mat, int &nrow, int &ncol)
{
	static const char *err = "The input should be a numeric matrix.";
	if (Rf_isMatrix(mat))
	{
		int *p = INTEGER(GET_DIM(mat));
		nrow = p[0]; ncol = p[1];
		if (TYPEOF(mat)!=REALSXP && !is_int(mat))
			throw_error_type(mat);
	} else if (is_sparse_seed(mat))
	{
		SEXP dm = GET_SLOT(mat, mkString("dim"));
		if (Rf_isNull(dm) || Rf_length(dm) != 2) Rf_error(err);
		int *p = INTEGER(dm);
		nrow = p[0]; ncol = p[1];
	} else
		throw_error_type(mat);
}



extern "C"
{

static inline double sq(double x) { return x*x; }

static inline double sum(const double x[], size_t n)
{
	double s = 0;
	for (size_t i=0; i < n; i++) s += x[i];
	return s;
}

static int block_idx_st = 0;

SEXP c_init_block()
{
	block_idx_st = 0;
	return R_NilValue;
}


// ====  a += b  ====

SEXP c_add_update(SEXP a, SEXP b)
{
	// 'a' should be a real vector or dense real matrix
	static const char *err_len = "a and b should have the same length.";
	// check
	if (TYPEOF(a) != REALSXP)
		Rf_error("The initial value should be numeric.");
	// add
	double *p = REAL(a);
	const size_t n = Rf_xlength(a);
	if (TYPEOF(b) == REALSXP)
	{
		if (Rf_xlength(a) != Rf_xlength(b)) Rf_error(err_len);
		const double *s = REAL(b);
		for (size_t i=0; i < n; i++) p[i] += s[i];
	} else if (is_int(b))
	{
		if (Rf_xlength(a) != Rf_xlength(b)) Rf_error(err_len);
		const int *s = INTEGER(b);
		for (size_t i=0; i < n; i++)
			if (s[i] != NA_INTEGER) p[i] += s[i]; else p[i] = NA_REAL;
	} else
		throw_error_type(b);
	// output
	return a;
}

// ====  rowSums & colSums  ====

SEXP c_rowSums(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pv = REAL(val);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) pv[j] += p[j];
		} else {
			FOR_LOOP_i_j  pv[j] += p[j];
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) pv[j] += p[j];
		} else {
			FOR_LOOP_i_j
				if (p[j] != NA_INTEGER) pv[j] += p[j]; else { pv[j] = NA_REAL; break; }
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i])) pv[M.nzi_r[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi  pv[M.nzi_r[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i];
					else pv[M.nzi_r[i]-1] = NA_REAL;
			}
		}
	}
	// output
	return val;
}

SEXP c_colSums(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Sum = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Sum);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s = 0;
				for (int j=0; j < nrow; j++) if (!ISNAN(p[j])) s += p[j];
				pv[i] = s;
			}
		} else {
			FOR_LOOP_i {
				double s = 0;
				for (int j=0; j < nrow; j++) s += p[j];
				pv[i] = s;
			}
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s = 0;
				for (int j=0; j < nrow; j++) if (p[j] != NA_INTEGER) s += p[j];
				pv[i] = s;
			}
		} else {
			FOR_LOOP_i {
				double s = 0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) s += p[j]; else { s = NA_REAL; break; }
				pv[i] = s;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		memset(pv, 0, sizeof(double)*ncol);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i])) pv[M.nzi_c[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi  pv[M.nzi_c[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER) pv[M.nzi_c[i]-1] += p[i];
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_c[i]-1] += p[i];
					else pv[M.nzi_c[i]-1] = NA_REAL;
			}
		}
	}
	// output
	UNPROTECT(1);
	return Sum;
}


// ====  rowProds & colProds  ====

SEXP c_rowProds(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	double *pv = REAL(val);
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) pv[j] *= p[j];
		} else {
			FOR_LOOP_i_j  pv[j] *= p[j];
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) pv[j] *= p[j];
		} else {
			FOR_LOOP_i_j
				if (p[j] != NA_INTEGER) pv[j] *= p[j]; else pv[j] = NA_REAL;
		}
	} else 
		throw_error_type(mat);
	// output
	return val;
}


SEXP c_colProds(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Prod = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Prod);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s = 1;
				for (int j=0; j < nrow; j++) if (!ISNAN(p[j])) s *= p[j];
				pv[i] = s;
			}
		} else {
			FOR_LOOP_i {
				double s = 1;
				for (int j=0; j < nrow; j++) s *= p[j];
				pv[i] = s;
			}
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s = 1;
				for (int j=0; j < nrow; j++) if (p[j] != NA_INTEGER) s *= p[j];
				pv[i] = s;
			}
		} else {
			FOR_LOOP_i {
				double s = 1;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) s *= p[j]; else { s = NA_REAL; break; }
				pv[i] = s;
			}
		}
	} else 
		throw_error_type(mat);
	// output
	return Prod;
}


// ====  rowMeans & colMeans  ====

SEXP c_rowMeans(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pv = REAL(val);
	int *pn = (int*)&pv[nrow];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) { pv[j] += p[j]; pn[j]++; }
		} else {
			FOR_LOOP_i_j  pv[j] += p[j];
			for (int j=0; j < nrow; j++) pn[j] += ncol;
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) { pv[j] += p[j]; pn[j]++; }
		} else {
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER) pv[j] += p[j];
				else { pv[j] = NA_REAL; break; }
			for (int j=0; j < nrow; j++) pn[j] += ncol;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		for (int j=0; j < nrow; j++) pn[j] += ncol;  // number of non-NaN
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1;
					if (!ISNAN(p[i])) pv[r] += p[i]; else pn[r]--;
				}
			} else {
				FOR_LOOP_Mi  pv[M.nzi_r[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER) pv[r] += p[i]; else pn[r]--;
				}
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER) pv[r] += p[i]; else pv[r] = NA_REAL;
				}
			}
		}
	}
	// output
	return val;
}

SEXP c_rowMeans_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 2;
	const double *sv = REAL(val);
	const int *sn = (int*)&sv[n];
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++) p[i] = sv[i] / sn[i];
	return ans;
}

SEXP c_colMeans(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Avg = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Avg);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0; int n=0;
				for (int j=0; j < nrow; j++)
					if (!ISNAN(p[j])) { s += p[j]; n++; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++) s += p[j];
				pv[i] = s / nrow;
			}
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0; int n=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { s += p[j]; n++; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) s += p[j]; else { s = NA_REAL; break; }
				pv[i] = s / nrow;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		memset(pv, 0, sizeof(double)*ncol);
		SEXP Num = PROTECT(NEW_INTEGER(ncol));
		int *pn = INTEGER(Num);
		for (int i=0; i < ncol; i++) pn[i] = nrow;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  if (!ISNAN(p[i]))
					pv[M.nzi_c[i]-1] += p[i]; else pn[M.nzi_c[i]-1]--;
			} else {
				FOR_LOOP_Mi  pv[M.nzi_c[i]-1] += p[i];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  if (p[i] != NA_INTEGER)
					pv[M.nzi_c[i]-1] += p[i]; else pn[M.nzi_c[i]-1]--;
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_c[i]-1] += p[i];
					else pv[M.nzi_c[i]-1] = NA_REAL;
			}
		}
		for (int i=0; i < ncol; i++) pv[i] /= pn[i];
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return Avg;
}


// ====  rowWeightedMeans & colWeightedMeans  ====

SEXP c_rowWMeans(SEXP mat, SEXP val, SEXP w, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const double *pw = REAL(w) + block_idx_st;
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pv = REAL(val), *pn = pv+nrow;
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j])) { pv[j] += p[j]*pw[i]; pn[j] += pw[i]; }
		} else {
			FOR_LOOP_i_j  pv[j] += p[j]*pw[i];
			const double w = sum(pw, ncol);
			for (int j=0; j < nrow; j++) pn[j] += w;
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
				{ pv[j] += p[j]*pw[i]; pn[j] += pw[i]; }
		} else {
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
				pv[j] += p[j]*pw[i]; else { pv[j] = NA_REAL; break; }
			const double w = sum(pw, ncol);
			for (int j=0; j < nrow; j++) pn[j] += w;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		const double w = sum(pw, ncol);
		for (int j=0; j < nrow; j++) pn[j] += w;
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i]-1, c = M.nzi_c[i]-1;
					if (!ISNAN(p[i])) pv[r] += p[i]*pw[c]; else pn[r] -= pw[c];
				}
			} else {
				FOR_LOOP_Mi  pv[M.nzi_r[i]-1] += p[i]*pw[M.nzi_c[i]-1];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i]-1, c = M.nzi_c[i]-1;
					if (p[i] != NA_INTEGER) pv[r] += p[i]*pw[c]; else pn[r] -= pw[c];
				}
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER) pv[M.nzi_r[i]-1] += p[i]*pw[M.nzi_c[i]-1];
					else pv[M.nzi_r[i]-1] = NA_REAL;
			}
		}
	}
	block_idx_st += ncol;
	// output
	return val;
}

SEXP c_rowWMeans_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 2;
	const double *sv = REAL(val), *sw = sv+n;
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++) p[i] = sv[i] / sw[i];
	return ans;
}

SEXP c_colWMeans(SEXP mat, SEXP w, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const double *pw = REAL(w);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Avg = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Avg), sum_w = 0;
	if (!na_rm)
		for (int j=0; j < nrow; j++) sum_w += pw[j];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, n=0;
				for (int j=0; j < nrow; j++)
					if (!ISNAN(p[j])) { s += p[j]*pw[j]; n += pw[j]; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++) s += p[j]*pw[j];
				pv[i] = s / sum_w;
			}
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, n=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { s += p[j]*pw[j]; n += pw[j]; }
				pv[i] = s / n;
			}
		} else {
			FOR_LOOP_i {
				double s=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) s += p[j]*pw[j]; else { s = NA_REAL; break; }
				pv[i] = s / sum_w;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		memset(pv, 0, sizeof(double)*ncol);
		SEXP Num = PROTECT(NEW_NUMERIC(ncol));
		double *pn = REAL(Num);
		for (int i=0; i < ncol; i++) pn[i] = sum_w;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i]-1, c = M.nzi_c[i]-1;
					if (!ISNAN(p[i])) pv[c] += p[i]*pw[r]; else pn[c] -= pw[r];
				}
			} else {
				FOR_LOOP_Mi  pv[M.nzi_c[i]-1] += p[i] * pw[M.nzi_r[i]-1];
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i]-1, c = M.nzi_c[i]-1;
					if (p[i] != NA_INTEGER) pv[c] += p[i]*pw[r]; else pn[c] -= pw[r];
				}
			} else {
				FOR_LOOP_Mi
					if (p[i] != NA_INTEGER)
						pv[M.nzi_c[i]-1] += p[i] * pw[M.nzi_r[i]-1];
					else pv[M.nzi_c[i]-1] = NA_REAL;
			}
		}
		for (int i=0; i < ncol; i++) pv[i] /= pn[i];
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return Avg;
}


// ====  rowVars & colVars  ====

SEXP c_rowVars(SEXP mat, SEXP val, SEXP narm, SEXP center)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const bool no_center = (Rf_isNull(center) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pS = REAL(val), *pS2 = pS+nrow;
	int *pN = (int*)&pS[2*nrow];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), v;
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i_j  if (!ISNAN(v=p[j]))
					{ pS[j] += v; pS2[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  { v = p[j]; pS[j] += v; pS2[j] += v*v; }
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i_j  if (!ISNAN(p[j]))
					{ v = p[j]-pC[j]; pS[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  { v = p[j]-pC[j]; pS[j] += v*v; }
			}
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat); double v;
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]; pS[j] += v; pS2[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]; pS[j] += v; pS2[j] += v*v; }
					else { pS[j] = pS2[j] = NA_REAL; break; }
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]-pC[j]; pS[j] += v*v; pN[j]++; }
			} else {
				FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
					{ v = p[j]-pC[j]; pS[j] += v*v; } else { pS[j] = NA_REAL; break; }
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		if (na_rm)
			for (int j=0; j < nrow; j++) pN[j] += ncol;
		if (!no_center)
		{
			const double *pC = REAL(center);
			for (int j=0; j < nrow; j++) pS[j] += ncol*pC[j]*pC[j];
		}
		int r;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata), v;
			if (no_center)
			{
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (!ISNAN(v=p[i])) { pS[r] += v; pS2[r] += v*v; } else pN[r]--;
					}
				} else {
					FOR_LOOP_Mi
						{ r = M.nzi_r[i]-1; v = p[i]; pS[r] += v; pS2[r] += v*v; }
				}
			} else {
				const double *pC = REAL(center);
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (!ISNAN(v=p[i])) pS[r] += v*v - 2*v*pC[r];
						else { pN[r]--; pS[r] -= pC[r]*pC[r]; }
					}
				} else {
					FOR_LOOP_Mi
						{ r = M.nzi_r[i]-1; v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
				}
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			double v;
			if (no_center)
			{
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v; pS2[r] += v*v; } else pN[r]--;
					}
				} else {
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pS[r] += v; pS2[r] += v*v; }
						else { pS[r] = pS2[r] = NA_REAL; }
					}
				}
			} else {
				const double *pC = REAL(center);
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
						else { pN[r]--; pS[r] -= pC[r]*pC[r]; }
					}
				} else {
					FOR_LOOP_Mi  {
						r = M.nzi_r[i] - 1;
						if (p[i] != NA_INTEGER) { v = p[i]; pS[r] += v*v - 2*v*pC[r]; }
						else pS[r] = NA_REAL;
					}
				}
			}
		}
	}
	// output
	if (!na_rm) for (int j=0; j < nrow; j++) pN[j] += ncol;
	return val;
}

SEXP c_rowVars_final(SEXP val, SEXP center)
{
	const bool no_center = Rf_isNull(center) == TRUE;
	const size_t n = Rf_xlength(val) / 3;
	const double *pS = REAL(val), *pS2 = pS+n;
	const int *pN = (const int*)&pS[2*n];
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	if (no_center)
	{
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? (pS2[i] - pS[i]*pS[i]/pN[i]) / (pN[i]-1) : NA_REAL;
	} else {
		for (size_t i=0; i < n; i++)
			p[i] = (pN[i]>1) ? pS[i] / (pN[i]-1) : NA_REAL;
	}
	return ans;
}

SEXP c_colVars(SEXP mat, SEXP narm, SEXP center)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const bool no_center = (Rf_isNull(center) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Var = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Var);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, s2=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (!ISNAN(v=p[j])) { s += v; s2 += v*v; n++; }
					pv[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, s2=0, v;
					for (int j=0; j < nrow; j++) { v = p[j]; s += v; s2 += v*v; }
					pv[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
				}
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (!ISNAN(p[j])) { v = p[j]-pC[i]; s += v*v; n++; }
					pv[i] = (n>1) ? s/(n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, v;
					for (int j=0; j < nrow; j++) { v = p[j]-pC[i]; s += v*v; }
					pv[i] = (nrow>1) ? s/(nrow-1) : NA_REAL;
				}
			}
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (no_center)
		{
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, s2=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; n++; }
					pv[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, s2=0, v;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; }
						else { s = s2 = NA_REAL; break; }
					pv[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
				}
			}
		} else {
			const double *pC = REAL(center);
			if (na_rm)
			{
				FOR_LOOP_i {
					double s=0, v; int n=0;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]-pC[i]; s += v*v; n++; }
					pv[i] = (n>1) ? s/(n-1) : NA_REAL;
				}
			} else {
				FOR_LOOP_i {
					double s=0, v;
					for (int j=0; j < nrow; j++)
						if (p[j] != NA_INTEGER) { v = p[j]-pC[i]; s += v*v; }
						else { s = NA_REAL; break; }
					pv[i] = (nrow>1) ? s/(nrow-1) : NA_REAL;
				}
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		int *pN=NULL;
		if (na_rm)
		{
			pN = INTEGER(PROTECT(NEW_INTEGER(ncol)));
			for (int i=0; i < ncol; i++) pN[i] = nrow;
		}
		// do sparse
		if (no_center)
		{
			memset(pv, 0, sizeof(double)*ncol);
			double *p2 = REAL(PROTECT(NEW_NUMERIC(ncol)));
			memset(p2, 0, sizeof(double)*ncol);
			if (TYPEOF(M.nzdata) == REALSXP)
			{
				double *p = REAL(M.nzdata), v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (!ISNAN(v=p[i])) { pv[c] += v; p2[c] += v*v; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						v = p[i]; pv[c] += v; p2[c] += v*v;
					}
				}
			} else {
				int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
				double v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v; p2[c] += v*v; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v; p2[c] += v*v; } else pv[c] = p2[c] = NA_REAL;
					}
				}
			}
			// final
			for (int i=0; i < ncol; i++)
			{
				int n = na_rm ? pN[i] : nrow;
				pv[i] = (n>1) ? (p2[i] - pv[i]*pv[i]/n) / (n-1) : NA_REAL;
			}
		} else {
			const double *pC = REAL(center);
			for (int i=0; i < ncol; i++) pv[i] = nrow*pC[i]*pC[i];
			if (TYPEOF(M.nzdata) == REALSXP)
			{
				double *p = REAL(M.nzdata), v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (!ISNAN(v=p[i])) pv[c] += v*v - 2*v*pC[c]; else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						v = p[i]; pv[c] += v*v - 2*v*pC[c];
					}
				}
			} else {
				int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
				double v;
				if (na_rm)
				{
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v*v - 2*v*pC[c]; } else pN[c]--;
					}
				} else {
					FOR_LOOP_Mi  {
						int c = M.nzi_c[i] - 1;
						if (p[i] != NA_INTEGER)
							{ v = p[i]; pv[c] += v*v - 2*v*pC[c]; } else pv[c] = NA_REAL;
					}
				}
			}
			// final
			for (int i=0; i < ncol; i++)
			{
				int n = na_rm ? pN[i] : nrow;
				pv[i] = (n>1) ? pv[i]/(n-1) : NA_REAL;
			}
		}
		// release
		UNPROTECT(no_center + na_rm);
	}
	// output
	UNPROTECT(1);
	return Var;
}


// ====  rowWeightedVars & colWeightedVars  ====

SEXP c_rowWVars(SEXP mat, SEXP val, SEXP w, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const double *pw = REAL(w) + block_idx_st;
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	double *pS = REAL(val), *pS2 = pS+nrow, *pSW = pS+2*nrow;
	int *pN = (int*)&pS[3*nrow];
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), v;
		if (na_rm)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j]))
			{
				pSW[j] += pw[i];
				v = p[j]*pw[i]; pS[j] += v; pS2[j] += v*v; pN[j]++;
			}
		} else {
			FOR_LOOP_i_j  {
				pSW[j] += pw[i];
				v = p[j]*pw[i]; pS[j] += v; pS2[j] += v*v;
			}
			for (int j=0; j < nrow; j++) pN[j] += ncol;
		}
	} else if (is_int(mat))
	{
		int *p = INTEGER(mat); double v;
		if (na_rm)
		{
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
			{
				pSW[j] += pw[i];
				v = p[j]*pw[i]; pS[j] += v; pS2[j] += v*v; pN[j]++;
			}
		} else {
			FOR_LOOP_i_j  if (p[j] != NA_INTEGER)
			{
				pSW[j] += pw[i];
				v = p[j]*pw[i]; pS[j] += v; pS2[j] += v*v;
			} else {
				pS[j] = pS2[j] = NA_REAL; break;
			}
			for (int j=0; j < nrow; j++) pN[j] += ncol;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		double sum_w = sum(pw, ncol);
		for (int j=0; j < nrow; j++) pSW[j] += sum_w;
		for (int j=0; j < nrow; j++) pN[j] += ncol;
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata), v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (!ISNAN(p[i]))
					{
						v = p[i]*pw[c]; pS[r] += v; pS2[r] += v*v;
					} else {
						pN[r]--; pSW[r] -= pw[c];
					}
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					v = p[i]*pw[c]; pS[r] += v; pS2[r] += v*v;
				}
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			double v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						v = p[i]*pw[c]; pS[r] += v; pS2[r] += v*v;
					} else {
						pN[r]--; pSW[r] -= pw[c];
					}
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						v = p[i]*pw[c]; pS[r] += v; pS2[r] += v*v;
					} else { pS[r] = pS2[r] = NA_REAL; }
				}
			}
		}
	}
	// output
	block_idx_st += ncol;
	return val;
}

SEXP c_rowWVars_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 4;
	const double *pS = REAL(val), *pS2 = pS+n, *pSW = pS+2*n;
	const int *pN = (const int*)&pS[3*n];
	SEXP ans = NEW_NUMERIC(n);
	double *p = REAL(ans);
	for (size_t i=0; i < n; i++)
	{
		double u = sq(pN[i] / pSW[i]);
		p[i] = (pN[i]>1) ? u*(pS2[i] - pS[i]*pS[i]/pN[i]) / (pN[i]-1) : NA_REAL;
	}
	return ans;
}

SEXP c_colWVars(SEXP mat, SEXP w, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	const double *pw = REAL(w);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP Var = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(Var);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, s2=0, sw=0, v; int n=0;
				for (int j=0; j < nrow; j++)
					if (!ISNAN(p[j]))
					{
						sw += pw[j];
						v = p[j]*pw[j]; s += v; s2 += v*v; n++;
					}
				pv[i] = (n>1) ? sq(n/sw)*(s2 - s*s/n) / (n-1) : NA_REAL;
			}
		} else {
			double u = sq(nrow / sum(pw, nrow));
			FOR_LOOP_i {
				double s=0, s2=0, v;
				for (int j=0; j < nrow; j++)
					{ v = p[j]*pw[j]; s += v; s2 += v*v; }
				pv[i] = (nrow>1) ? u*(s2 - s*s/nrow) / (nrow-1) : NA_REAL;
			}
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, s2=0, sw=0, v; int n=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER)
					{
						sw += pw[j];
						v = p[j]*pw[j]; s += v; s2 += v*v; n++;
					}
				pv[i] = (n>1) ? sq(n/sw)*(s2 - s*s/n) / (n-1) : NA_REAL;
			}
		} else {
			double u = sq(nrow / sum(pw, nrow));
			FOR_LOOP_i {
				double s=0, s2=0, v;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { v = p[j]*pw[j]; s += v; s2 += v*v; }
					else { s = s2 = NA_REAL; break; }
				pv[i] = (nrow>1) ? u*(s2 - s*s/nrow) / (nrow-1) : NA_REAL;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		int *pN=NULL;
		if (na_rm)
		{
			pN = INTEGER(PROTECT(NEW_INTEGER(ncol)));
			for (int i=0; i < ncol; i++) pN[i] = nrow;
		}
		// do sparse
		memset(pv, 0, sizeof(double)*ncol);
		double *p2 = REAL(PROTECT(NEW_NUMERIC(ncol*2)));
		memset(p2, 0, sizeof(double)*ncol);
		double *sw = p2 + ncol;
		double sum_w = sum(pw, nrow);
		for (int i=0; i < ncol; i++) sw[i] = sum_w;

		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata), v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (!ISNAN(p[i]))
					{
						v = p[i]*pw[r]; pv[c] += v; p2[c] += v*v;
					} else {
						pN[c]--; sw[c] -= pw[r];
					}
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					v = p[i]*pw[r]; pv[c] += v; p2[c] += v*v;
				}
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			double v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						v = p[i]*pw[r]; pv[c] += v; p2[c] += v*v;
					} else {
						pN[c]--; sw[c] -= pw[r];
					}
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1, r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						v = p[i]*pw[r]; pv[c] += v; p2[c] += v*v;
					} else {
						pv[c] = p2[c] = NA_REAL;
					}
				}
			}
		}
		// final
		for (int i=0; i < ncol; i++)
		{
			int n = na_rm ? pN[i] : nrow;
			double u = sq(n / sw[i]);
			pv[i] = (n>1) ? u*(p2[i] - pv[i]*pv[i]/n) / (n-1) : NA_REAL;
		}
		// release
		UNPROTECT(1 + (int)na_rm);
	}
	// output
	UNPROTECT(1);
	return Var;
}


// ====  scRowMeanVar & scColMeanVar  ====

// also: SEXP c_rowVars(SEXP mat, SEXP val, SEXP narm, SEXP center)

SEXP c_rowMeanVar_final(SEXP val)
{
	const size_t n = Rf_xlength(val) / 3;
	const double *pS = REAL(val), *pS2 = pS+n;
	const int *pN = (const int*)&pS[2*n];
	SEXP ans = Rf_allocMatrix(REALSXP, n, 2);
	double *pM = REAL(ans), *pV = pM + n;
	for (size_t i=0; i < n; i++)
	{
		const int n = pN[i];
		pM[i] = (n > 0) ? pS[i]/n : NA_REAL;
		pV[i] = (n > 1) ? (pS2[i] - pS[i]*pS[i]/n) / (n-1) : NA_REAL;
	}
	return ans;
}

SEXP c_colMeanVar(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// output variable
	SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, ncol, 2));
	double *pv = REAL(ans), *pv2 = pv + ncol;
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, s2=0, v; int n=0;
				for (int j=0; j < nrow; j++)
					if (!ISNAN(v=p[j])) { s += v; s2 += v*v; n++; }
				pv[i]  = (n>0) ? s/n : NA_REAL;
				pv2[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
			}
		} else {
			FOR_LOOP_i {
				double s=0, s2=0, v;
				for (int j=0; j < nrow; j++) { v = p[j]; s += v; s2 += v*v; }
				pv[i]  = s / nrow;
				pv2[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
			}
		}
	} else if (is_int(mat))
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double s=0, s2=0, v; int n=0;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; n++; }
				pv[i]  = (n>0) ? s/n : NA_REAL;
				pv2[i] = (n>1) ? (s2 - s*s/n) / (n-1) : NA_REAL;
			}
		} else {
			FOR_LOOP_i {
				double s=0, s2=0, v;
				for (int j=0; j < nrow; j++)
					if (p[j] != NA_INTEGER) { v = p[j]; s += v; s2 += v*v; }
					else { s = s2 = NA_REAL; break; }
				pv[i]  = s / nrow;
				pv2[i] = (nrow>1) ? (s2 - s*s/nrow) / (nrow-1) : NA_REAL;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		// initialize
		memset(pv, 0, sizeof(double)*ncol*2);
		int *pN=NULL;
		if (na_rm)
		{
			pN = INTEGER(PROTECT(NEW_INTEGER(ncol)));
			for (int i=0; i < ncol; i++) pN[i] = nrow;
		}
		// do sparse
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			double *p = REAL(M.nzdata), v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					if (!ISNAN(v=p[i])) { pv[c] += v; pv2[c] += v*v; } else pN[c]--;
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					v = p[i]; pv[c] += v; pv2[c] += v*v;
				}
			}
		} else {
			int *p = INTEGER(M.nzdata);  // INTSXP or LGLSXP
			double v;
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					if (p[i] != NA_INTEGER)
						{ v = p[i]; pv[c] += v; pv2[c] += v*v; } else pN[c]--;
				}
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					if (p[i] != NA_INTEGER)
						{ v = p[i]; pv[c] += v; pv2[c] += v*v; } else pv[c] = pv2[c] = NA_REAL;
				}
			}
		}
		// final
		for (int i=0; i < ncol; i++)
		{
			int n = na_rm ? pN[i] : nrow;
			double sum = pv[i];
			pv[i]  = sum / n;
			pv2[i] = (n>1) ? (pv2[i] - sum*sum/n) / (n-1) : NA_REAL;
		}
		// release
		UNPROTECT((int)na_rm);
	}
	// output
	UNPROTECT(1);
	return ans;
}

}
