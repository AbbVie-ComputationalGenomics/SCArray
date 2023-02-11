// ===============================================================
//
// SCArray R package
// Copyright (C) 2023    Xiuwen Zheng
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


extern "C"
{
// ====  rowMins & colMins  ====

SEXP c_rowMins(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	double *pv = REAL(val);
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j {
			double v = p[j], &d = pv[j];
			if (!ISNAN(v))
			{
				if (!ISNAN(d) && (v < d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j {
			int v = p[j]; double &d = pv[j];
			if (v != NA_INTEGER)
			{
				if (!ISNAN(d) && (v < d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(nrow)));
		memset(pn, 0, sizeof(int)*nrow);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double v = p[i], &d = pv[r];
				if (!ISNAN(v))
				{
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double &d = pv[r];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < nrow; i++)
			if (pn[i]<ncol && !ISNAN(pv[i]) && 0<pv[i]) pv[i] = 0;
		UNPROTECT(1);
	}
	// output
	return val;
}

SEXP c_colMins(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	SEXP ans = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(ans);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i {
			double d=R_PosInf, v;
			for (int j=0; j < nrow; j++)
			{
				if (!ISNAN(v=p[j]))
				{
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) { d = NA_REAL; break; }
				}			
			}
			pv[i] = d;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i {
			double d=R_PosInf; int v;
			for (int j=0; j < nrow; j++)
			{
				if ((v=p[j]) != NA_INTEGER)
				{
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) { d = NA_REAL; break; }
				}			
			}
			pv[i] = d;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(ncol)));
		memset(pn, 0, sizeof(int)*ncol);
		for (int i=0; i < ncol; i++) pv[i] = R_PosInf;
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double v = p[i], &d = pv[c];
				if (!ISNAN(v))
				{
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double &d = pv[c];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d) && (v < d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < ncol; i++)
			if (pn[i]<nrow && !ISNAN(pv[i]) && 0<pv[i]) pv[i] = 0;
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return ans;
}


// ====  rowMaxs & colMaxs  ====

SEXP c_rowMaxs(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	double *pv = REAL(val);
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j {
			double v = p[j], &d = pv[j];
			if (!ISNAN(v))
			{
				if (!ISNAN(d) && (v > d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j {
			int v = p[j]; double &d = pv[j];
			if (v != NA_INTEGER)
			{
				if (!ISNAN(d) && (v > d)) d = v;
			} else {
				if (!na_rm) d = NA_REAL;
			}			
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(nrow)));
		memset(pn, 0, sizeof(int)*nrow);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double v = p[i], &d = pv[r];
				if (!ISNAN(v))
				{
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double &d = pv[r];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < nrow; i++)
			if (pn[i]<ncol && !ISNAN(pv[i]) && pv[i]<0) pv[i] = 0;
		UNPROTECT(1);
	}
	// output
	return val;
}

SEXP c_colMaxs(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	// get nrow and ncol
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	SEXP ans = PROTECT(NEW_NUMERIC(ncol));
	double *pv = REAL(ans);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i {
			double d=R_NegInf, v;
			for (int j=0; j < nrow; j++)
			{
				if (!ISNAN(v=p[j]))
				{
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) { d = NA_REAL; break; }
				}			
			}
			pv[i] = d;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i {
			double d=R_NegInf; int v;
			for (int j=0; j < nrow; j++)
			{
				if ((v=p[j]) != NA_INTEGER)
				{
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) { d = NA_REAL; break; }
				}			
			}
			pv[i] = d;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(ncol)));
		memset(pn, 0, sizeof(int)*ncol);
		for (int i=0; i < ncol; i++) pv[i] = R_NegInf;
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double v = p[i], &d = pv[c];
				if (!ISNAN(v))
				{
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double &d = pv[c];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d) && (v > d)) d = v;
				} else {
					if (!na_rm) d = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < ncol; i++)
			if (pn[i]<nrow && !ISNAN(pv[i]) && pv[i]<0) pv[i] = 0;
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return ans;
}


// ====  rowRanges & colRanges  ====

SEXP c_rowRanges(SEXP mat, SEXP val, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	// get nrow and ncol
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	double *pv1 = REAL(val), *pv2 = pv1+nrow;
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j {
			double v=p[j], &d1=pv1[j], &d2=pv2[j];
			if (!ISNAN(v))
			{
				if (!ISNAN(d1) && (v < d1)) d1 = v;
				if (!ISNAN(d2) && (v > d2)) d2 = v;
			} else {
				if (!na_rm) d1 = d2 = NA_REAL;
			}			
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j {
			int v = p[j]; double &d1=pv1[j], &d2=pv2[j];
			if (v != NA_INTEGER)
			{
				if (!ISNAN(d1) && (v < d1)) d1 = v;
				if (!ISNAN(d2) && (v > d2)) d2 = v;
			} else {
				if (!na_rm) d1 = d2 = NA_REAL;
			}			
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(nrow)));
		memset(pn, 0, sizeof(int)*nrow);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double v = p[i], &d1 = pv1[r], &d2 = pv2[r];
				if (!ISNAN(v))
				{
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) d1 = d2 = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int r = M.nzi_r[i] - 1;
				pn[r] ++;
				double &d1 = pv1[r], &d2 = pv2[r];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) d1 = d2 = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < nrow; i++)
		{
			if (pn[i] < ncol)
			{
				if (!ISNAN(pv1[i]) && 0<pv1[i]) pv1[i] = 0;
				if (!ISNAN(pv2[i]) && 0>pv2[i]) pv2[i] = 0;
			}
		}
		UNPROTECT(1);
	}
	// output
	return val;
}

SEXP c_colRanges(SEXP mat, SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, ncol, 2));
	double *pv1 = REAL(ans), *pv2 = pv1+ncol;
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i {
			double d1=R_PosInf, d2=R_NegInf, v;
			for (int j=0; j < nrow; j++)
			{
				if (!ISNAN(v=p[j]))
				{
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) { d1 = d2 = NA_REAL; break; }
				}			
			}
			pv1[i] = d1; pv2[i] = d2;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i {
			double d1=R_PosInf, d2=R_NegInf; int v;
			for (int j=0; j < nrow; j++)
			{
				if ((v=p[j]) != NA_INTEGER)
				{
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) { d1 = d2 = NA_REAL; break; }
				}			
			}
			pv1[i] = d1; pv2[i] = d2;
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		int *pn = INTEGER(PROTECT(NEW_INTEGER(ncol)));
		memset(pn, 0, sizeof(int)*ncol);
		for (int i=0; i < ncol; i++)
			pv1[i] = R_PosInf, pv2[i] = R_NegInf;
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double v = p[i], &d1 = pv1[c], &d2 = pv2[c];
				if (!ISNAN(v))
				{
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) d1 = d2 = NA_REAL;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi {
				int c = M.nzi_c[i] - 1;
				pn[c] ++;
				double &d1 = pv1[c], &d2 = pv2[c];
				if (p[i] != NA_INTEGER)
				{
					double v = p[i];
					if (!ISNAN(d1) && (v < d1)) d1 = v;
					if (!ISNAN(d2) && (v > d2)) d2 = v;
				} else {
					if (!na_rm) d1 = d2 = NA_REAL;
				}
			}
		}
		// check zero
		for (int i=0; i < ncol; i++)
		{
			if (pn[i] < nrow)
			{
				if (!ISNAN(pv1[i]) && 0<pv1[i]) pv1[i] = 0;
				if (!ISNAN(pv2[i]) && 0>pv2[i]) pv2[i] = 0;
			}
		}
		UNPROTECT(1);
	}
	// output
	UNPROTECT(1);
	return ans;
}

}
