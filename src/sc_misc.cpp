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
#include <cstdint>
#include <vector>


extern "C"
{

// ====  scNumSplit  ====

inline static SEXP range(int64_t st, int64_t ed)
{
	const int64_t i32max = 2147483647;
	SEXP ans;
	if ((st <= i32max) && (ed <= i32max))
	{
		ans = NEW_INTEGER(2);
		int *p = INTEGER(ans);
		p[0] = st; p[1] = ed;
	} else {
		ans = NEW_NUMERIC(2);
		double *p = REAL(ans);
		p[0] = st; p[1] = ed;
	}
	return ans;
}

SEXP c_get_split(SEXP totnum, SEXP numworker, SEXP buffer)
{
	// check totnum
	const double ntot_d = Rf_asReal(totnum);
	if (!R_FINITE(ntot_d) || (ntot_d < 0))
		Rf_error("the total number must be >= 0.");
	const int64_t ntot = (int64_t)ntot_d;
	if (ntot <= 0) return NEW_LIST(0);
	// check numworker
	const int nw = Rf_asInteger(numworker);
	if (nw < 0) Rf_error("# of workers must be > 0.");
	if (nw <= 1)
	{
		SEXP ans = PROTECT(NEW_LIST(1));
		SET_ELEMENT(ans, 0, range(1, ntot));
		UNPROTECT(1);
		return ans;
	}
	// split
	int64_t *pbuf = (int64_t*)REAL(buffer);
	const double r = (double)ntot / nw;
	double st = 0;
	for (int i=0; i < nw; i++, st+=r) pbuf[i] = (int64_t)round(st);
	pbuf[nw] = ntot;
	// output
	int n = 0;
	for (int i=0; i < nw; i++)
		if (pbuf[i] < pbuf[i+1]) n++;
	SEXP ans = PROTECT(NEW_LIST(n));
	n = 0;
	for (int i=0; i < nw; i++)
	{
		if (pbuf[i] < pbuf[i+1])
			SET_ELEMENT(ans, n++, range(pbuf[i]+1, pbuf[i+1]));
	}
	UNPROTECT(1);
	return ans;
}


// ====  scRowAutoGrid & scColAutoGrid  ====

inline static int max_0(int v) { return (v >= 0) ? v : 0; }

SEXP c_sparse_blocksize(SEXP tot_nbyte, SEXP range_max, SEXP nnzero,
	SEXP buffer)
{
	double nb_f = Rf_asReal(tot_nbyte);
	if (!R_FINITE(nb_f)) nb_f = 1;
	int64_t nbyte = (int64_t)nb_f;
	if (nbyte <= 0) nbyte = 1;
	// adjust memory to # of non-zero
	nbyte /= (4 + 4 + 8);  // icol, irow, value
	nbyte /= 2;            // need buffer
	if (nbyte <= 0) nbyte = 1;
	// max. of block size
	double rg_max = Rf_asReal(range_max);
	if (!R_FINITE(rg_max)) rg_max = 1;
	int n_max = (int)rg_max;
	if (n_max < 0)
		n_max = 2147483647;
	else if (n_max == 0)
		n_max = 1;
	// calculate
	const int num = Rf_length(nnzero);
	int *p_nz = INTEGER(nnzero), *p_buf = INTEGER(buffer);
	int n_len = 0;
	for (int i=0; i < num; )
	{
		int old_i = i;
		int64_t m = 0;
		for (; (i < num) && ((m + max_0(p_nz[i])) <= nbyte); i++)
			m += max_0(p_nz[i]);
		if (old_i == i) i++;
		int n = i - old_i;
		if (n > n_max) { i -= n - n_max; n = n_max; }
		p_buf[n_len++] = n;
	}
	// ouptut
	SEXP ans = PROTECT(NEW_INTEGER(n_len));
	int *p = INTEGER(ans);
	for (int i=0; i < n_len; i++) p[i] = p_buf[i];
	UNPROTECT(1);
	return ans;
}


// ====  rowsum & colsum  ====

SEXP c_row_sum_grp(SEXP mat, SEXP ii, SEXP n_grp, SEXP narm)
{
	const int n_group = Rf_asInteger(n_grp);
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	const int *p_ii = INTEGER(ii);
	// output variable
	SEXP ans = PROTECT(Rf_allocMatrix(REALSXP, n_group, ncol));
	double *pv = REAL(ans);
	memset(pv, 0, sizeof(double)*n_group*ncol);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				for (int j=0; j < nrow; j++)
					if (ISNAN(p[j])) pv[p_ii[j]] += p[j];
				pv += n_group;
			}
		} else {
			FOR_LOOP_i {
				for (int j=0; j < nrow; j++) pv[p_ii[j]] += p[j];
				pv += n_group;
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				for (int j=0; j < nrow; j++)
					if (p[j] == NA_INTEGER) pv[p_ii[j]] += p[j];
				pv += n_group;
			}
		} else {
			FOR_LOOP_i {
				for (int j=0; j < nrow; j++)
					if (p[j] == NA_INTEGER) pv[p_ii[j]] += p[j];
					else { pv[p_ii[j]] = NA_REAL; break; }
				pv += n_group;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					if (ISNAN(p[i])) pv[p_ii[r] + c*size_t(n_group)] += p[i];
				}
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					pv[p_ii[r] + c*size_t(n_group)] += p[i];
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					if (p[i] != NA_INTEGER) pv[p_ii[r] + c*size_t(n_group)] += p[i];
				}
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					size_t k = p_ii[r] + c*size_t(n_group);
					if (p[i] != NA_INTEGER) pv[k] += p[i]; else pv[k] = NA_REAL;
				}
			}
		}
	}
	// output
	UNPROTECT(1);
	return ans;
}

SEXP c_col_sum_grp(SEXP mat, SEXP out_mat, SEXP ii, SEXP start, SEXP n_grp,
	SEXP narm)
{
	const bool na_rm = (Rf_asLogical(narm) == TRUE);
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	const int st = INTEGER(start)[1] - 1;
	const int *p_ii = INTEGER(ii) + st;
	// output variable
	double *pv = REAL(out_mat);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double *pv2 = pv + p_ii[i] * size_t(nrow);
				for (int j=0; j < nrow; j++)
					if (ISNAN(p[j])) pv2[j] += p[j];
			}
		} else {
			FOR_LOOP_i {
				double *pv2 = pv + p_ii[i] * size_t(nrow);
				for (int j=0; j < nrow; j++) pv2[j] += p[j];
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_rm)
		{
			FOR_LOOP_i {
				double *pv2 = pv + p_ii[i] * size_t(nrow);
				for (int j=0; j < nrow; j++)
					if (p[j] == NA_INTEGER) pv2[j] += p[j];
			}
		} else {
			FOR_LOOP_i {
				double *pv2 = pv + p_ii[i] * size_t(nrow);
				for (int j=0; j < nrow; j++)
					if (p[j] == NA_INTEGER) pv2[j] += p[j]; else pv2[j] = NA_REAL;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					if (ISNAN(p[i])) pv[p_ii[c] * size_t(nrow) + r] += p[i];
				}
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					pv[p_ii[c] * size_t(nrow) + r] += p[i];
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_rm)
			{
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					if (p[i] != NA_INTEGER) pv[p_ii[c] * size_t(nrow) + r] += p[i];
				}
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1, c = M.nzi_c[i] - 1;
					size_t k = p_ii[c] * size_t(nrow) + r;
					if (p[i] != NA_INTEGER) pv[k] += p[i]; else pv[k] = NA_REAL;
				}
			}
		}
	}
	// output
	return out_mat;
}


// ====  row_nnzero & col_nnzero  ====

SEXP c_row_nnzero(SEXP mat, SEXP val, SEXP na)
{
	const int na_counted = LOGICAL(na)[0];
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	int *pv = INTEGER(val);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_counted == TRUE)
		{
			FOR_LOOP_i_j  if (ISNAN(p[j]) || (p[j]!=0)) pv[j]++;
		} else if (na_counted == FALSE)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j]) && (p[j]!=0)) pv[j]++;
		} else {
			FOR_LOOP_i_j  {
				if (!ISNAN(p[j]))
				{
					if ((p[j] != 0) && (pv[j] != NA_INTEGER)) pv[j]++;
				} else
					pv[j] = NA_INTEGER;
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_counted == TRUE)
		{
			FOR_LOOP_i_j  if ((p[j]==NA_INTEGER) || (p[j]!=0)) pv[j]++;
		} else if (na_counted == FALSE)
		{
			FOR_LOOP_i_j  if ((p[j]!=NA_INTEGER) && (p[j]!=0)) pv[j]++;
		} else {
			FOR_LOOP_i_j  {
				if (p[j] != NA_INTEGER)
				{
					if ((p[j] != 0) && (pv[j] != NA_INTEGER)) pv[j]++;
				} else
					pv[j] = NA_INTEGER;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			if (na_counted == TRUE)
			{
				FOR_LOOP_Mi
					if (ISNAN(p[i]) || (p[i]!=0)) pv[M.nzi_r[i]-1]++;
			} else if (na_counted == FALSE)
			{
				FOR_LOOP_Mi
					if (!ISNAN(p[i]) && (p[i]!=0)) pv[M.nzi_r[i]-1]++;
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1;
					if (!ISNAN(p[i]))
					{
						if ((p[i] != 0) && (pv[r] != NA_INTEGER)) pv[r]++;
					} else
						pv[r] = NA_INTEGER;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_counted == TRUE)
			{
				FOR_LOOP_Mi
					if ((p[i]==NA_INTEGER) || (p[i]!=0)) pv[M.nzi_r[i]-1]++;
			} else if (na_counted == FALSE)
			{
				FOR_LOOP_Mi
					if ((p[i]!=NA_INTEGER) && (p[i]!=0)) pv[M.nzi_r[i]-1]++;
			} else {
				FOR_LOOP_Mi  {
					int r = M.nzi_r[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						if ((p[i] != 0) && (pv[r] != NA_INTEGER)) pv[r]++;
					} else
						pv[r] = NA_INTEGER;
				}
			}
		}
	}
	// output
	return val;
}

SEXP c_col_nnzero(SEXP mat, SEXP na)
{
	const int na_counted = LOGICAL(na)[0];
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	SEXP ans = PROTECT(NEW_INTEGER(ncol));
	int *pv = INTEGER(ans);
	memset(pv, 0, sizeof(int)*ncol);  // set all FALSE
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		if (na_counted == TRUE)
		{
			FOR_LOOP_i_j  if (ISNAN(p[j]) || (p[j]!=0)) pv[i]++;
		} else if (na_counted == FALSE)
		{
			FOR_LOOP_i_j  if (!ISNAN(p[j]) && (p[j]!=0)) pv[i]++;
		} else {
			FOR_LOOP_i_j  {
				if (!ISNAN(p[j]))
				{
					if ((p[j] != 0) && (pv[i] != NA_INTEGER)) pv[i]++;
				} else
					pv[i] = NA_INTEGER;
			}
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		if (na_counted == TRUE)
		{
			FOR_LOOP_i_j  if ((p[j]==NA_INTEGER) || (p[j]!=0)) pv[i]++;
		} else if (na_counted == FALSE)
		{
			FOR_LOOP_i_j  if ((p[j]!=NA_INTEGER) && (p[j]!=0)) pv[i]++;
		} else {
			FOR_LOOP_i_j  {
				if (p[j] != NA_INTEGER)
				{
					if ((p[j] != 0) && (pv[i] != NA_INTEGER)) pv[i]++;
				} else
					pv[i] = NA_INTEGER;
			}
		}
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			if (na_counted == TRUE)
			{
				FOR_LOOP_Mi
					if (ISNAN(p[i]) || (p[i]!=0)) pv[M.nzi_c[i]-1]++;
			} else if (na_counted == FALSE)
			{
				FOR_LOOP_Mi
					if (!ISNAN(p[i]) && (p[i]!=0)) pv[M.nzi_c[i]-1]++;
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					if (!ISNAN(p[i]))
					{
						if ((p[i] != 0) && (pv[c] != NA_INTEGER)) pv[c]++;
					} else
						pv[c] = NA_INTEGER;
				}
			}
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			if (na_counted == TRUE)
			{
				FOR_LOOP_Mi
					if ((p[i]==NA_INTEGER) || (p[i]!=0)) pv[M.nzi_c[i]-1]++;
			} else if (na_counted == FALSE)
			{
				FOR_LOOP_Mi
					if ((p[i]!=NA_INTEGER) && (p[i]!=0)) pv[M.nzi_c[i]-1]++;
			} else {
				FOR_LOOP_Mi  {
					int c = M.nzi_c[i] - 1;
					if (p[i] != NA_INTEGER)
					{
						if ((p[i] != 0) && (pv[c] != NA_INTEGER)) pv[c]++;
					} else
						pv[c] = NA_INTEGER;
				}
			}
		}
	}
	// output
	UNPROTECT(1);
	return ans;
}


// ====  rowAnyNAs & colAnyNAs  ====

SEXP c_row_anyNA(SEXP mat, SEXP flag)
{
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	int *pv = LOGICAL(flag);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j  if (ISNAN(p[j])) pv[j] = TRUE;
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j  if (p[j] == NA_INTEGER) pv[j] = TRUE;
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi
				if (ISNAN(p[i])) pv[M.nzi_r[i] - 1] = TRUE;
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi
				if (p[i] == NA_INTEGER) pv[M.nzi_r[i] - 1] = TRUE;
		}
	}
	// output
	return flag;
}

SEXP c_col_anyNA(SEXP mat)
{
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	SEXP ans = PROTECT(NEW_LOGICAL(ncol));
	int *pv = LOGICAL(ans);
	memset(pv, 0, sizeof(int)*ncol);  // set all FALSE
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		const double *p = REAL(mat);
		FOR_LOOP_i_j  if (ISNAN(p[j])) { pv[i] = TRUE; break; }
	} else if (TYPEOF(mat) == INTSXP)
	{
		const int *p = INTEGER(mat);
		FOR_LOOP_i_j  if (p[j] == NA_INTEGER) { pv[i] = TRUE; break; }
	} else if (is_sparse_seed(mat))
	{
		SparseMatrix M(mat);
		if (TYPEOF(M.nzdata) == REALSXP)
		{
			const double *p = REAL(M.nzdata);
			FOR_LOOP_Mi
				if (ISNAN(p[i])) pv[M.nzi_c[i] - 1] = TRUE;
		} else {
			const int *p = INTEGER(M.nzdata);  // INTSXP
			FOR_LOOP_Mi
				if (p[i] == NA_INTEGER) pv[M.nzi_c[i] - 1] = TRUE;
		}
	}
	// output
	UNPROTECT(1);
	return ans;
}


// ====  rowCollapse & colCollapse  ====

static std::vector< std::vector<int> > row_map;
static int idx_col = 0;

SEXP c_rowCollapse_init(SEXP idx, SEXP dim)
{
	const int nrow = INTEGER(dim)[0];
	const int ncol = INTEGER(dim)[1];
	const int *pIdx = INTEGER(idx), nIdx = Rf_length(idx);
	row_map.clear();
	row_map.resize(ncol);
	int ii = 0;
	for (int i=0; i < nrow; i++)
	{
		int j = pIdx[ii++];
		if (ii >= nIdx) ii = 0;
		if (0 < j && j <= ncol)
			row_map[j-1].push_back(i);
	}
	idx_col = 0;
	return R_NilValue;
}

SEXP c_rowCollapse_done()
{
	row_map.clear();
	return R_NilValue;
}

SEXP c_rowCollapse(SEXP mat, SEXP val)
{
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), *pv = REAL(val);
		FOR_LOOP_i {
			const std::vector<int> &r = row_map[idx_col + i];
			std::vector<int>::const_iterator it = r.begin();
			for (; it != r.end(); it++)
				pv[*it] = p[*it];
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		int *p = INTEGER(mat), *pv = INTEGER(val);
		FOR_LOOP_i {
			const std::vector<int> &r = row_map[idx_col + i];
			std::vector<int>::const_iterator it = r.begin();
			for (; it != r.end(); it++)
				pv[*it] = p[*it];
		}
	}
	idx_col += ncol;
	// output
	return val;
}

// ====

static int idx_n = 0, idx_i = 0;
static int *idx_p = NULL;

SEXP c_colCollapse_init(SEXP idx)
{
	idx_n = Rf_length(idx);
	idx_i = idx_col = 0;
	idx_p = (idx_n > 0) ? INTEGER(idx) : &NA_INTEGER;
	return R_NilValue;
}

SEXP c_colCollapse(SEXP mat, SEXP val)
{
	int nrow, ncol;
	get_mat_size(mat, nrow, ncol);
	// do
	if (TYPEOF(mat) == REALSXP)
	{
		double *p = REAL(mat), *pv = REAL(val);
		FOR_LOOP_i {
			int j = idx_p[idx_i++];
			if (idx_i >= idx_n) idx_i = 0;
			pv[idx_col+i] = (j != NA_INTEGER) ? p[j-1] : NA_REAL;
		}
	} else if (TYPEOF(mat) == INTSXP)
	{
		int *p = INTEGER(mat), *pv = INTEGER(val);
		FOR_LOOP_i {
			int j = idx_p[idx_i++];
			if (idx_i >= idx_n) idx_i = 0;
			pv[idx_col+i] = (j != NA_INTEGER) ? p[j-1] : j;
		}
	}
	idx_col += ncol;
	// output
	return val;
}

}
