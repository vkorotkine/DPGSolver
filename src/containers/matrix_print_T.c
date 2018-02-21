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

#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_tol.h"

// Static function declarations ************************************************************************************* //

///\{ \name The number of decimal places to display.
#define N_DEC 4
///\}

/** \brief Check whether either of the matrix extents is zero.
 *  \return `true` if yes; `false` otherwise. */
static bool check_Matrix_extents_zero_T
	(const ptrdiff_t ext_0, ///< First extent.
	 const ptrdiff_t ext_1  ///< Second extent.
	);

// Interface functions ********************************************************************************************** //

void print_Matrix_T_tol (const struct Matrix_T*const a, const Real tol)
{
	MAYBE_UNUSED(tol);
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	if (check_Matrix_extents_zero_T(ext_0,ext_1))
		printf("Called print_Matrix_T for input with extents: [%td,%td].\n\n",ext_0,ext_1);

	const Type* data = a->data;

	char format_d[10];
	char format_i[10];
	sprintf(format_d,"%s%de%s","% .",N_DEC," ");
	sprintf(format_i,"%s%dd%s","% ",N_DEC+7," ");

	switch (a->layout) {
	case 'R':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const Type val = *data++;
#ifdef TYPE_RC
	#if TYPE_RC == TYPE_REAL
				if (isnan(val) || (fabs(val) > tol))
					printf(format_d,val);
				else
					printf(format_i,0);
	#elif TYPE_RC == TYPE_COMPLEX
				print_real(val);
				print_imag(val);
	#endif
#else
				printf(format_i,val);
#endif
			}
			printf("\n");
		}
		printf("\n");
		break;
	case 'C':
		for (ptrdiff_t i = 0; i < ext_0; i++) {
			for (ptrdiff_t j = 0; j < ext_1; j++) {
				const Type val = data[i+ext_0*j];
#ifdef TYPE_RC
	#if TYPE_RC == TYPE_REAL
				if (isnan(val) || (fabs(val) > tol))
					printf(format_d,val);
				else
					printf(format_i,0);
	#elif TYPE_RC == TYPE_COMPLEX
				print_real(val);
				print_imag(val);
	#endif
#else
				printf(format_i,val);
#endif
			}
			printf("\n");
		}
		printf("\n");
		break;
	default:
		EXIT_ERROR("Unsupported: %c",a->layout);
		break;
	}
}

void print_const_Matrix_T_tol (const struct const_Matrix_T*const a, const Real tol)
{
	print_Matrix_T_tol((struct Matrix_T*)a,tol);
}

void print_Matrix_T (const struct Matrix_T*const a)
{
	print_Matrix_T_tol(a,EPS);
}

void print_const_Matrix_T (const struct const_Matrix_T*const a)
{
	print_Matrix_T((const struct Matrix_T*)a);
}

void print_to_file_const_Matrix_T (FILE*const file, const struct const_Matrix_T*const a)
{
	const ptrdiff_t ext_0 = a->ext_0,
	                ext_1 = a->ext_1;

	const Type* data = a->data;

	static const char*const format_d = " % .4e"; MAYBE_UNUSED(format_d);
	static const char*const format_i = " % 11d";  MAYBE_UNUSED(format_i);
	switch (a->layout) {
	case 'R':
		for (int i = 0; i < ext_0; ++i) {
			for (int j = 0; j < ext_1; ++j) {
				const Type val = *data++;
#ifdef TYPE_RC
	#if TYPE_RC == TYPE_REAL
				if (isnan(val) || (fabs(val) > 0.0))
					fprintf(file,format_d,val);
				else
					fprintf(file,format_i,0);
	#elif TYPE_RC == TYPE_COMPLEX
EXIT_ADD_SUPPORT; UNUSED(file); UNUSED(val);
	#endif
#else
				fprintf(file,format_i,val);
#endif
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
		break;
	case 'C':
		for (int i = 0; i < ext_0; ++i) {
			for (int j = 0; j < ext_1; ++j) {
				const Type val = data[i+ext_0*j];
#ifdef TYPE_RC
	#if TYPE_RC == TYPE_REAL
				if (isnan(val) || (fabs(val) > 0.0))
					fprintf(file,format_d,val);
				else
					fprintf(file,format_i,0);
	#elif TYPE_RC == TYPE_COMPLEX
EXIT_ADD_SUPPORT; UNUSED(file); UNUSED(val);
	#endif
#else
				fprintf(file,format_i,val);
#endif
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
		break;
	default:
		EXIT_ERROR("Unsupported: %c",a->layout);
		break;
	}
}

void print_to_file_Matrix_T (FILE*const file, const struct Matrix_T*const a)
{
	print_to_file_const_Matrix_T(file,(struct const_Matrix_T*)a);
}

#if TYPE_RC == TYPE_COMPLEX
void print_real (const double complex val)
{
	static char format_d[10],
	            format_i[10];

	static bool format_set = false;
	if (!format_set) {
		format_set = true;
		sprintf(format_d,"%s%de","% .",N_DEC);
		sprintf(format_i,"%s%dd","% ",N_DEC+7);
	}

	static const double tol = EPS;
	const double val_r = creal(val);
	if (isnan(val_r) || (fabs(val_r) > tol))
		printf(format_d,val_r);
	else
		printf(format_i,0);
}

void print_imag (const double complex val)
{
	static char format_dp[10],
	            format_dm[10],
	            format_ip[10];

	static bool format_set = false;
	if (!format_set) {
		format_set = true;
		sprintf(format_dp,"%s%de%s","+%.",N_DEC,"i ");
		sprintf(format_dm,"%s%de%s","% .",N_DEC,"i ");
		sprintf(format_ip,"%s%dd%s","+% ",N_DEC+6,"i ");
	}

	static const double tol_i = EPS*EPS*EPS;
	const double val_i = cimag(val);
	if (isnan(val_i) || (fabs(val_i) > tol_i)) {
		if (val_i >= 0.0)
			printf(format_dp,val_i);
		else
			printf(format_dm,val_i);
	} else {
		printf(format_ip,0);
	}
}
#endif

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static bool check_Matrix_extents_zero_T (const ptrdiff_t ext_0, const ptrdiff_t ext_1)
{
	if (ext_0 == 0 || ext_1 == 0)
		return true;
	return false;
}
