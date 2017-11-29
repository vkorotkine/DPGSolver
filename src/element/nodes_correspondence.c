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

#include "nodes_correspondence.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "element.h"
#include "nodes.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for a \ref const_Vector_i holding the face correspondence indices for the input permutation
 *         index.
 *  \return See brief. */
static const struct const_Vector_i* constructor_face_corr
	(const struct const_Nodes* nodes, ///< \ref const_Nodes.
	 const int ind_perm,              ///< The index of the permutation.
	 const int e_type                 ///< \ref Element::type.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_Vector_i* constructor_nodes_face_corr
	(const int d, const int p, const int node_type, const int s_type)
{
	assert(d >= 0 && d <= 2); // faces only.
	constructor_Nodes_fptr constructor_Nodes = get_constructor_Nodes_by_super_type(s_type);

	const struct const_Nodes* nodes = constructor_Nodes(d,p,node_type); // destructed

	const ptrdiff_t n_perm = get_n_perm_corr(d,s_type);

	struct Multiarray_Vector_i* face_corr = constructor_empty_Multiarray_Vector_i(false,1,&n_perm); // returned

	const int e_type = compute_elem_from_super_type(s_type,d);
	for (int i = 0; i < n_perm; ++i)
		face_corr->data[i] = (struct Vector_i*) constructor_face_corr(nodes,i,e_type);

	destructor_const_Nodes(nodes);

	return (const struct const_Multiarray_Vector_i*) face_corr;
}

ptrdiff_t get_n_perm_corr (const int d, const int s_type)
{
	switch (d) {
	case 0:
		return POINT_N_PERM;
		break;
	case 1:
		return LINE_N_PERM;
		break;
	case 2:
		switch (s_type) {
			case ST_TP: return QUAD_N_PERM; break;
			case ST_SI: return TRI_N_PERM;  break;
			default: EXIT_ERROR("Unsupported: %d.\n",s_type); break;
		}
	default:
		EXIT_ERROR("Unsupported: %d.\n",d);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Reverse the entries of the input data.
static void reverse_entries
	(const int n, ///< The number of entries.
	 int* data    ///< The data.
	);

/// \brief Swap the two input data blocks.
static void swap_blocks
	(const int n, ///< The number of entries in both of the blocks.
	 int* data_0, ///< The data of block 0.
	 int* data_1  ///< The data of block 1.
	);

static const struct const_Vector_i* constructor_face_corr
	(const struct const_Nodes* nodes, const int ind_perm, const int e_type)
{
	/// \todo Needs clean-up.
	const ptrdiff_t d  = nodes->rst->ext_1;
	const ptrdiff_t Nn = nodes->rst->ext_0;
	const double* rst  = nodes->rst->data;

	struct Vector_i* face_corr = constructor_empty_Vector_i(Nn);
	int* fc_data = face_corr->data;

	int i, jMax, iInd;

	switch (d) {
	case 0:
		fc_data[0] = 0;
		break;
	case 1:
		for (int i = 0; i < Nn; i++)
			fc_data[i] = i;
		switch (ind_perm) {
		case 0:
			; // Do nothing
			break;
		case 1:
			reverse_entries((int)Nn,&fc_data[0]);
			break;
		default:
			EXIT_ERROR("Unsupported: %d.\n",ind_perm);
			break;
		}
		break;
	case 2:
		if (e_type == QUAD) {
			int sqrtNn = (int)sqrt(Nn);

			switch(ind_perm) {
			case 0: case 1: case 2: case 3:
				for (int i = 0; i < Nn; i++)
					fc_data[i] = i;
				break;
			case 4: case 5: case 6: case 7:
				for (int i = 0; i < sqrtNn; i++) {
					iInd = i*sqrtNn;
					for (int j = 0; j < sqrtNn; j++)
						fc_data[iInd+j] = i+j*sqrtNn;
				}
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_perm);
				break;
			}

			// Reverse entries of 1D-blocks if necessary
			switch (ind_perm) {
			case 0: case 2: case 4: case 5:
				; // Do nothing
				break;
			case 1: case 3: case 6: case 7:
				for (int i = 0; i < sqrtNn; ++i)
					reverse_entries(sqrtNn,&fc_data[i*sqrtNn]);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_perm);
				break;
			}

			// Swap 1D-blocks if necessary
			switch (ind_perm) {
			case 0: case 1: case 4: case 6:
				; // Do nothing
				break;
			case 2: case 3: case 5: case 7:
				for (int i = 0; i < sqrtNn/2; ++i)
					swap_blocks(sqrtNn,&fc_data[i*sqrtNn],&fc_data[(sqrtNn-1-i)*sqrtNn]);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",ind_perm);
				break;
			}
		} else if (e_type == TRI) {
			assert(nodes->has_symms == true);

			const ptrdiff_t Ns = nodes->s->ext_0;
			const int* symms = nodes->s->data;

			int j, k, kMax, iInd, subOrder[3], fc_data_swap3[3], fc_data_swap[Nn], Foundn[Nn], IndX[Nn];
			double       DY[Nn*Nn];

			for (i = 0; i < Nn; i++)
				fc_data[i] = i;

			// Swap entries if necessary
			switch(ind_perm) {
				default: // default cases 0, 1, 2
					; // Do nothing
					break;
				case 3: case 4: case 5:
					for (i = 0; i < Nn; i++) {
						iInd = (int)(i*Nn);
						for (j = 0; j < Nn; j++) {
							DY[iInd+j] = fabs(rst[Nn+i]-rst[Nn+j]);
						}
					}

					for (i = 0; i < Nn; i++)
						Foundn[i] = 0;

					for (i = 0; i < Nn; i++) {
						if (!Foundn[i]) {
							iInd = (int)(i*Nn);
							kMax = 0;
							for (j = 0; j < Nn; j++) {
								if (!Foundn[j] && i != j && DY[iInd+j] < 10*EPS)
									IndX[kMax++] = j;
							}
							for (k = 0; k < kMax; k++) {
								if (fabs(rst[i]+rst[IndX[k]]) < 1e3*EPS) {
									Foundn[i] = 1;
									Foundn[IndX[k]] = 1;
									fc_data_swap[i] = IndX[k];
									fc_data_swap[IndX[k]] = i;

									break;
								}
							}
							if (kMax == 0) {
								Foundn[i] = 1;
								fc_data_swap[i] = i;
							}
						}
					}

					for (i = 0; i < Nn; i++) {
						if (Foundn[i] == 0)
							printf("Error: Did not find all nodes in get_face_ordering (TRI).\n"), exit(1);
					}

					for (i = 0; i < Nn; i++)
						fc_data[i] = fc_data_swap[i];
					break;
			}

			// Rotate entries of 3-symmetry blocks if necessary
			switch(ind_perm) {
				case 0:
				case 5:
					subOrder[0] = 0; subOrder[1] = 1; subOrder[2] = 2;
					break;
				case 1:
				case 3:
					subOrder[0] = 1; subOrder[1] = 2; subOrder[2] = 0;
					break;
				case 2:
				case 4:
					subOrder[0] = 2; subOrder[1] = 0; subOrder[2] = 1;
					break;
				default:
					EXIT_ERROR("Unsupported: %d\n",ind_perm);
					break;
			}

			iInd = 0;
			for (i = 0; i < Ns; i++) {
				if (i) iInd += symms[i-1];
				jMax = symms[i];
				if (jMax == 3) {
					for (j = 0; j < jMax; j++)
						fc_data_swap3[j] = fc_data[iInd+subOrder[j]];
					for (j = 0; j < jMax; j++)
						fc_data[iInd+j] = fc_data_swap3[j];
				}
				// Setting the 1-symmetry orbit (if present) is redundant as the node position remains unchanged.
			}
		} else {
			printf("Error: Unsupported e_type in 3D in get_face_ordering.\n"), exit(1);
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %td.\n",d);
		break;
	}

	return (const struct const_Vector_i*) face_corr;
}

// Level 1 ********************************************************************************************************** //

static void reverse_entries (const int n, int* data)
{
	for (int i = 0; i < n/2; ++i) {
		const int tmp = data[i];
		data[i] = data[n-1-i];
		data[n-1-i] = tmp;
	}
}

static void swap_blocks (const int n, int* data_0, int* data_1)
{
	for (int i = 0; i < n; ++i) {
		const int tmp = data_0[i];
		data_0[i] = data_1[i];
		data_1[i] = tmp;
	}
}
