/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */

#include "test_support_math_functions.h"

#include <assert.h>

#include "macros.h"
#include "definitions_tol.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void z_yxpz_dcc (const int n, const double* x, const double complex* y, double complex* z)
{
	for (int i = 0; i < n; ++i)
		z[i] += y[i]*x[i];
}

double norm_diff_petsc_Mat (Mat A0, Mat A1)
{
	PetscInt m[2] = { 0, 0, },
	         n[2] = { 0, 0, };
	MatGetSize(A0,&m[0],&n[0]);
	MatGetSize(A1,&m[1],&n[1]);

	assert(m[0] == m[1]);
	assert(n[0] == n[1]);

	double norm_num = 0.0,
	       norm_den = 0.0;

	assert(m[0] == n[0]);
	const int i_max = m[0];
	for (int i = 0; i < i_max; ++i) {
		int n_cols[2] = { 0, 0, };
		const PetscInt* cols[2]    = { NULL, NULL, };
		const PetscScalar* vals[2] = { NULL, NULL, };

		MatGetRow(A0,i,&n_cols[0],&cols[0],&vals[0]);
		MatGetRow(A1,i,&n_cols[1],&cols[1],&vals[1]);

		if (n_cols[0] != n_cols[1])
			EXIT_ERROR("Differing number of entries in row %d (%d, %d).\n",i,n_cols[0],n_cols[1]);
		const double diff = norm_diff_d(n_cols[0],vals[0],vals[1],"Inf");
		if (diff > norm_num)
			norm_num = diff;

		const double norm_row[2] = { norm_d(n_cols[0],vals[0],"Inf"), norm_d(n_cols[1],vals[1],"Inf"), };
		const double max = ( norm_row[0] > norm_row[1] ? norm_row[0] : norm_row[1] );
		if (max > norm_den)
			norm_den = max;
		assert(max < 1e30);
printf("%d %e %e\n",i,norm_num,norm_den);

		MatRestoreRow(A0,i,&n_cols[0],&cols[0],&vals[0]);
		MatRestoreRow(A1,i,&n_cols[1],&cols[1],&vals[1]);
	}

	return ( fabs(norm_den) > 1e2*EPS ? norm_num/norm_den : norm_num );
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
