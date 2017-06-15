// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "memory_constructors_matrix.h"

#include <stdlib.h>

#include "Macros.h"
#include "Parameters.h"

#include "matrix_structs.h"
#include "setup_operators_support.h"
#include "element_functions.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for allocation of matrix structs.
 *
 *	Comments:
 *		The final allocation of the 1-level dereferenced structs is performed by the appropriate setup_operators
 *		function.
 *
 *	Notation:
 *
 *	References:
 */

struct S_MATRIX *constructor1_mat_D (const char layout, const size_t NRows, const size_t NCols, const size_t NColsSub)
{
	struct S_MATRIX *A = calloc(1,sizeof *A); // returned

	A->layout = layout;
	A->format = 'D';
	A->NRows  = NRows;
	A->NCols  = NCols;
	A->NColsSub = NColsSub;
	A->values = calloc(NRows*NCols , sizeof *(A->values)); // keep

	return A;
}

struct S_MATRIX **constructor2_mat (size_t const N0)
{
	struct S_MATRIX **A = calloc(N0 , sizeof *A);
	return A;
}

struct S_MATRIX ***constructor3_mat (size_t const N0, size_t const N1)
{
	struct S_MATRIX ***A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor2_mat(N1);
	return A;
}

struct S_MATRIX ****constructor4_mat (size_t const N0, size_t const N1, size_t const N2)
{
	struct S_MATRIX ****A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor3_mat(N1,N2);
	return A;
}

struct S_MATRIX *****constructor5_mat (size_t const N0, size_t const N1, size_t const N2, size_t const N3)
{
	struct S_MATRIX *****A = calloc(N0 , sizeof *A);
	for (size_t i = 0; i < N0; i++)
		A[i] = constructor4_mat(N1,N2,N3);
	return A;
}


void constructor_move2_mat (char const layout, const char format, const size_t NRows, const size_t NCols,
                            double *const values, struct S_MATRIX **A)
{
	/*
	 *	Comments:
	 *		Assumed that a single matrix is returned.
	 */

	struct S_MATRIX *A1 = calloc(1,sizeof *A1); // returned

	A1->layout = layout;
	A1->format = format;
	A1->NRows  = NRows;
	A1->NCols  = NCols;
	A1->values = values;

	A[0] = A1;
}

void constructor_move4_mat (char const layout, char const format,
                            unsigned int const *const NRows1, unsigned int const *const NCols1,
							unsigned int const *const *const NRows2, unsigned int const *const *const NCols2,
							double ****A_d, struct S_MATRIX ****A, struct S_OP_RANGE *const op_range)
{
	if (format != 'D')
		EXIT_UNSUPPORTED; // Add support

	if (op_range->type_op == 'V') {
		set_operator_ranges(op_range,'S');
		for (size_t PS = op_range->PSMin; PS <= op_range->PSMax; PS++) {
			op_range->PS = PS;
			set_operator_ranges(op_range,'I');
			set_operator_ranges(op_range,'V');
			for (size_t Pb = op_range->PbMin; Pb <= op_range->PbMax; Pb++) {
			for (size_t vh = op_range->vhMin; vh < op_range->vhMax; vh++) {
				struct S_MATRIX *A1 = calloc(1,sizeof *A1); // keep

				A1->layout = layout;
				A1->format = format;
				A1->NRows  = NRows1[Pb];
				A1->NCols  = NCols1[PS];
				A1->values = A_d[PS][Pb][vh];

				A[PS][Pb][vh] = A1;
			}}
		}
	} else if (op_range->type_op == 'F') {
		set_operator_ranges(op_range,'S');
		for (size_t PS = op_range->PSMin; PS <= op_range->PSMax; PS++) {
			op_range->PS = PS;
			set_operator_ranges(op_range,'I');
			for (size_t Pb = op_range->PbMin; Pb <= op_range->PbMax; Pb++) {
			for (size_t f  = 0; f < op_range->ELEMENT->Nf; f++) {
				size_t const Eclass   = get_Eclass(op_range->EType);
				size_t const IndFType = get_IndFType(Eclass,f);

				op_range->f = f;
				set_operator_ranges(op_range,'F');
				for (size_t fh = op_range->fhMin; fh < op_range->fhMax; fh++) {
					size_t const ft = f*NFREFMAX+fh;

					struct S_MATRIX *A1 = calloc(1,sizeof *A1); // keep

					A1->layout = layout;
					A1->format = format;

					switch (op_range->e_to_e) {
					case RvCf_op_r:
						A1->NRows = NRows1[Pb];
						A1->NCols = NCols2[PS][IndFType];
						break;
					case RfCv_op_r:
						A1->NRows = NRows2[Pb][IndFType];
						A1->NCols = NCols1[PS];
						break;
					case RfCf_op_r:
						A1->NRows = NRows2[Pb][IndFType];
						A1->NCols = NCols2[PS][IndFType];
						break;
					case RvCv_op_r:
						// fallthrough (Invalid for face operators)
					default:
						EXIT_UNSUPPORTED;
						break;
					}
					A1->values = A_d[PS][Pb][ft];

					A[PS][Pb][ft] = A1;
				}
			}}
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void constructor_move5_mat (char const layout, char const format, unsigned int const *const NRows,
                            unsigned int const *const NCols, double *****A_d, struct S_MATRIX *****A,
                            struct S_OP_RANGE *const op_range)
{
	if (format != 'D')
		EXIT_UNSUPPORTED; // Add support

	if (op_range->type_op == 'V') {
		set_operator_ranges(op_range,'S');
		for (size_t PS = op_range->PSMin; PS <= op_range->PSMax; PS++) {
			op_range->PS = PS;
			set_operator_ranges(op_range,'I');
			set_operator_ranges(op_range,'V');
			set_operator_ranges(op_range,'D');
			for (size_t Pb = op_range->PbMin; Pb <= op_range->PbMax; Pb++) {
			for (size_t vh = op_range->vhMin; vh < op_range->vhMax; vh++) {
			for (size_t d  = op_range->dMin;  d  <= op_range->dMax;  d++)  {
				struct S_MATRIX *A1 = calloc(1,sizeof *A1); // keep

				A1->layout = layout;
				A1->format = format;
				A1->NRows  = NRows[Pb];
				A1->NCols  = NCols[PS];
				A1->values = A_d[PS][Pb][vh][d];

				A[PS][Pb][vh][d] = A1;
			}}}
		}
	} else if (op_range->type_op == 'F') {
		EXIT_UNSUPPORTED; // Add support
	} else {
		EXIT_UNSUPPORTED;
	}
}
