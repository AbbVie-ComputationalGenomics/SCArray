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
#include <vector>


extern "C"
{
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
