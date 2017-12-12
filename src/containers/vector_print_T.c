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
/// \file

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"

#include "matrix_print_T.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void print_Vector_T_tol (const struct Vector_T*const a, const Real tol)
{
	MAYBE_UNUSED(tol);
	const ptrdiff_t ext = a->ext_0;

	const Type* data = a->data;

	for (ptrdiff_t i = 0; i < ext; i++) {
#ifdef TYPE_RC
		const Type val = *data++;
	#if TYPE_RC == TYPE_REAL
		printf("% .3e ",( (isnan(val) || (fabs(val) > tol)) ? val : 0.0 ));
	#elif TYPE_RC == TYPE_COMPLEX
		print_real(val);
		print_imag(val);
	#endif
#else
		printf("% 12d ",*data++);
#endif
		if (!((i+1)%8))
			printf(" ...\n");
	}
	printf("\n\n");
}

void print_const_Vector_T_tol (const struct const_Vector_T*const a, const Real tol)
{
	print_Vector_T_tol((struct Vector_T*)a,tol);
}

void print_Vector_T (const struct Vector_T*const a)
{
	print_Vector_T_tol(a,EPS);
}

void print_const_Vector_T (const struct const_Vector_T*const a)
{
	print_Vector_T((const struct Vector_T*)a);
}

#ifndef TYPE_RC
void fprint_const_Vector_T (FILE* file, const int n_tab, const struct const_Vector_T* a)
{
	for (int i = 0; i < n_tab; ++i)
		fprintf(file,"\t");

	const ptrdiff_t ext_0 = a->ext_0;
	for (ptrdiff_t i = 0; i < ext_0; ++i)
		fprintf(file," % 3d",a->data[i]);

	fprintf(file,"\n");
}

void fprint_Vector_T (FILE* file, const int n_tab, struct Vector_T* a)
{
	fprint_const_Vector_T(file,n_tab,(const struct const_Vector_T*)a);
}
#endif

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
