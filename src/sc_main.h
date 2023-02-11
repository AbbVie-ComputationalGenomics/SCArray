// ===============================================================
//
// SCArray R package
// Copyright (C) 2022-2023    Xiuwen Zheng
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

#ifndef H_SCARRAY
#define H_SCARRAY


#include <R.h>
#include <Rdefines.h>


// Define C Macro

#define FOR_LOOP_i    \
	for (int i=0; i < ncol; i++, p+=nrow)
#define FOR_LOOP_i_j    \
	for (int i=0; i < ncol; i++, p+=nrow) for (int j=0; j < nrow; j++)
#define FOR_LOOP_Mi    for (int i=0; i < M.nnzero; i++)


/// SparseArraySeed
struct SparseMatrix
{
	int nnzero;
	const int *nzi_r, *nzi_c;
	SEXP nzdata;
	/// constructor	
	SparseMatrix(SEXP mat);
};


// Return true if mat is a SparseArraySeed object
bool is_sparse_seed(SEXP mat);
// Get # of rows and columns
void get_mat_size(SEXP mat, int &nrow, int &ncol);


#endif    // H_SCARRAY
