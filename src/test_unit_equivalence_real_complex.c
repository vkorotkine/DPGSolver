// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_equivalence_real_complex.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "fluxes_inviscid.h"
#include "fluxes_inviscid_c.h"

//#include "array_print.h"

/*
 *	Purpose:
 *		Test equivalence of output between real and complex versions of the same functions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static unsigned int compare_flux_inviscid(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                          const unsigned int Neq, double *Wr)
{
	unsigned int pass = 0;

	unsigned int   i, iMax, NnTotal, Nvar;
	double         *Fr, *Fctr;
	double complex *Wc, *Fc;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	Wc   = malloc(NnTotal*Nvar  * sizeof *Wc);   // free
	Fr   = malloc(NnTotal*d*Neq * sizeof *Fr);   // free
	Fc   = malloc(NnTotal*d*Neq * sizeof *Fc);   // free
	Fctr = malloc(NnTotal*d*Neq * sizeof *Fctr); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		Wc[i] = Wr[i];

	flux_inviscid(Nn,Nel,Wr,Fr,d,Neq);
	flux_inviscid_c(Nn,Nel,Wc,Fc,d,Neq);

	for (i = 0, iMax = NnTotal*d*Neq; i < iMax; i++)
		Fctr[i] = creal(Fc[i]);

	pass = 0;
	if (array_norm_diff_d(NnTotal*d*Neq,Fr,Fctr,"Inf") < EPS)
		pass = 1, TestDB.Npass++;
	
	free(Wc);
	free(Fr);
	free(Fc);
	free(Fctr);

	return pass;
}

void test_unit_equivalence_real_complex(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *		Real functions input.
	 *
	 *	Expected Output:
	 *		No difference between outputs.
	 *
	 */

	unsigned int Nn, Nel, d, Neq;
	double       *Wr;

	for (d = 1; d <= 3; d++) {
		Neq  = d+2;

		Wr   = initialize_W(&Nn,&Nel,d); // free
		pass = compare_flux_inviscid(Nn,Nel,d,Neq,Wr);

		if (d == 1) {
			//     0         10        20        30        40        50
			printf("equivalence_flux_inviscid (d = %d):               ",d);
		} else {
			printf("                          (d = %d):               ",d);
		}
		test_print(pass);

		free(Wr);
	}

	// flux_ROE
	// flux_LF
	// explicit_VOLUME_info
}
