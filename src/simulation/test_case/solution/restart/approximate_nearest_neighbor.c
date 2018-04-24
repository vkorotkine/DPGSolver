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
 *
 *  In preparation for a possible extension to larger integer types, we write many of the `int`eger data types using the
 *  generic "Index" type. This may however require replacement with templated container names.
 *
 *  The implementation is based on the minimalist [ANN algorithm of Chan]. However, as this algorithm operates only on
 *  integer types, the scalar values are first converted to their
 *
 *  <!-- References: -->
 *  [ANN algorithm of Chan]: ann/Chan2006_A_Minimalist's_Implementation_of_an_ANN_Algorithm_in_Fixed_Dimensions.pdf
 */

#include "approximate_nearest_neighbor.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "macros.h"
#include "definitions_core.h"

#include "matrix.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

#define SINGLE ///< Definition to use `int` type for functions below.

/// \{ \name Templating-related definitions.
#ifdef SINGLE
	#define Index     int
	#define INDEX_MIN 0
	#define INDEX_MAX (INT_MAX/4)
#else
	#define Index     long long
	#define INDEX_MIN 0
	#define INDEX_MAX (LLONG_MAX/4)
#endif
#define Real double
#define REAL_MAX DBL_MAX
/// \}

Index shift; ///< The shift to use. Global in this file as it must be used by \todo ref cmp_shuffle passed to qsort.

/** The accuracy tolerance for the approximate nearest neighbor search. Chosen based on the limitations of the
 *  conversion to `int`. */
#define EPS_ANN 1e-9

#define BASE_SEED 31415926535 ///< Random number base seed.

#define POW2_R(x) (((Real) (x))*((Real) (x))) ///< 'R'eal pow(x,2).

/// \brief Container for a node represented in binary.
struct Node {
	Index xyz[DIM]; ///< The coordinates with type `Index`.
	int index;      ///< The index.
};

/// \brief Container for 'S'hift-'S'huffle-'S'ort approximate nearest neighbor information.
struct SSS {
	ptrdiff_t ext_b; ///< The number of entries in \ref SSS::b.
	ptrdiff_t ext_s; ///< The number of entries in \ref SSS::s.

	const struct Node* b;  ///< The background nodes.
	const struct Node* s;  ///< The search nodes.
};

/// \brief Container for the converging node list and current search node.
struct SSS_c {
	ptrdiff_t n_b;  ///< The number of background nodes.
	struct Node* b; ///< The background nodes.
	struct Node* s; ///< The pointer to the search node.

	struct Node l; ///< Lower bounding node.
	struct Node u; ///< Upper bounding node.

	// ANN-related parameters
	Real r2;          ///< The minimum squared Euclidian distance from the search node to computed background nodes.
	int ind_ann;      ///< The current index of the approximate nearest neighbor.
	struct Node* ann; ///< The pointer to approximate nearest neighbor node.
};

/** \brief Constructor for a \ref SSS container.
 *  \return See brief. */
static const struct SSS* constructor_SSS
	(const struct Input_ANN*const ann_i ///< Standard.
	);

/// \brief Destructor for a \ref SSS container.
static void destructor_SSS
	(const struct SSS*const sss ///< Standard.
	);

/** \brief Perform the binary search for the input node in the list of background nodes.
 *
 *  ANN-related parameters of \ref SSS_c are set based on the current best guess before returning.
 */
static void SSS_query
	(struct SSS_c*const sss ///< Standard.
	);

// Interface functions ********************************************************************************************** //

const struct const_Vector_i* constructor_ann_indices (const struct Input_ANN*const ann_i)
{
	const ptrdiff_t n_b = ann_i->nodes_b->ext_0,
	                n_s = ann_i->nodes_s->ext_0;

	srand48(BASE_SEED+n_b+n_s+DIM); // seed the random number generator.
	shift = (Index) (drand48()*INDEX_MAX);

	const struct SSS*const sss = constructor_SSS(ann_i); // destructed

	struct Vector_i* ann_indices = constructor_empty_Vector_i(n_s); // returned
	for (int n = 0; n < n_s; ++n) {
		struct SSS_c sss_c =
			{ .n_b = n_b,
			  .b   = (struct Node*) sss->b,
			  .s   = (struct Node*) &(sss->s[n]),
			  .r2  = REAL_MAX,
			};
		SSS_query(&sss_c);
		ann_indices->data[n] = sss_c.ann->index;
	}
	destructor_SSS(sss);

	return (struct const_Vector_i*) ann_indices;
}

void destructor_Input_ANN (struct Input_ANN*const ann_info)
{
	destructor_const_Matrix_d(ann_info->nodes_b);
	destructor_const_Matrix_d(ann_info->nodes_s);
	free(ann_info);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Print an array of \ref Node\*s.
static void print_nodes
	(const ptrdiff_t n,             ///< Array length.
	 const struct Node*const nodes, ///< The array.
	 const char*const name          ///< Name to display.
	);

/** \brief Constructor for a Matrix holding the values of the minimum and maximum node coordinates in each of the xyz
 *         coordinate directions; the last two rows are used to hold the average ((min+max)/2) and the difference
 *         (max-min).
 *  \return See brief. */
static const struct const_Matrix_d* constructor_xyz_min_max
	(const struct Input_ANN*const ann_i ///< Standard.
	);

/** \brief Shuffle order comparison function for nodes.
 *  \return Standard required return for comparator function of qsort.
 *
 *  [Chan] (p.4) proved: p <~ q iff x_j <= y_j, j is the index of the 'm'ost 's'ignificant 'b'it of (x xor y) for j
 *  \f$\in\f$ d. The <~ operator, called the shuffle order, is used to denote "less than" comparison of the shuffle of
 *  two points written in binary. See [Chan] section 2.
 *
 *  <!-- References: -->
 *  [Chan]: ann/Chan2006_A_Minimalist's_Implementation_of_an_ANN_Algorithm_in_Fixed_Dimensions.pdf
 */
static int cmp_shuffle
	(Index* p, ///< The 1st node (represented in binary).
	 Index* q  ///< The 2nd node (represented in binary).
	);

/** \brief Compute the distance between the background node of index n and the search node and update
 *         ANN-related members of \ref SSS_c if the current node is closer than all previous. */
static void compute_distance_and_update
	(const ptrdiff_t n,     ///< The index of \ref SSS_c::b.
	 struct SSS_c*const sss ///< Standard.
	);

/** \brief Return the squared distance from the search node to the current bounding box.
 *  \return See brief. */
Real compute_r2_to_box
	(const struct SSS_c*const sss ///< Standard.
	);

static const struct SSS* constructor_SSS (const struct Input_ANN*const ann_i)
{
	enum { display = 0, };

	struct SSS*const sss = calloc(1,sizeof *sss); // free

	const ptrdiff_t n_b = ann_i->nodes_b->ext_0,
	                n_s = ann_i->nodes_s->ext_0;

	sss->ext_b = n_b;
	sss->ext_s = n_s;
	struct Node*const b = calloc((size_t)n_b,sizeof *(sss->b)); // moved
	struct Node*const s = calloc((size_t)n_s,sizeof *(sss->s)); // moved

	const struct const_Matrix_d*const xyz_mm = constructor_xyz_min_max(ann_i); // destructed
	const double*const xyz_avg = get_row_const_Matrix_d(2,xyz_mm),
	            *const xyz_h   = get_row_const_Matrix_d(3,xyz_mm);

	const double* data_b = ann_i->nodes_b->data,
	            * data_s = ann_i->nodes_s->data;
	for (int n = 0; n < n_b; ++n) {
		b[n].index = n;
		for (int d = 0; d < DIM; ++d) {
			const double rst = 2.0/xyz_h[d]*(*data_b-xyz_avg[d]);
			b[n].xyz[d] = (Index) (0.5*((1.0-rst)*INDEX_MIN+(1.0+rst)*INDEX_MAX));
			++data_b;
		}
	}
	for (int n = 0; n < n_s; ++n) {
		for (int d = 0; d < DIM; ++d) {
			const double rst = 2.0/xyz_h[d]*(*data_s-xyz_avg[d]);
			s[n].xyz[d] = (Index) (0.5*((1.0-rst)*INDEX_MIN+(1.0+rst)*INDEX_MAX));
			++data_s;
		}
	}
	destructor_const_Matrix_d(xyz_mm);

	if (display)
		print_nodes((int)n_b,b,"background");
	qsort((void*)b,(size_t)n_b,sizeof *b,(int (*)(const void *, const void *))cmp_shuffle);
	if (display) {
		print_nodes((int)n_b,b,"background - sorted");
		print_nodes((int)n_s,s,"search");
	}

	sss->b = b; // free
	sss->s = s; // free

	return sss;
}

static void destructor_SSS (const struct SSS*const sss)
{
	free((void*)sss->b);
	free((void*)sss->s);
	free((void*)sss);
}

static void SSS_query (struct SSS_c*const sss)
{
	const ptrdiff_t n_b = sss->n_b;
	if (n_b == 0)
		return;

	compute_distance_and_update(n_b/2,sss);
	const Real r2 = sss->r2;
//UNUSED(r2);

	if (n_b == 1 || compute_r2_to_box(sss)*POW2_R(1+EPS_ANN) > r2)
//	if (n_b == 1)
		return;

	struct Node*const b = sss->b,
	           *const s = sss->s;

	if (cmp_shuffle(s->xyz,b[n_b/2].xyz) < 0) { // p.3 line 4
		sss->n_b = n_b/2; // Chan p.3 line 5 (binary search in lower half)
		SSS_query(sss);

		if (cmp_shuffle(sss->u.xyz,b[n_b/2].xyz) > 0) { // Chan p.3 line 6
			sss->b   = &b[n_b/2+1];
			sss->n_b = n_b-(n_b/2+1);
			SSS_query(sss);
		}
	} else {
		sss->b   = &b[n_b/2+1];   // p.3 line 8 (binary search in upper half)
		sss->n_b = n_b-(n_b/2+1);
		SSS_query(sss);

		if (cmp_shuffle(sss->l.xyz,b[n_b/2].xyz) < 0) { // Chan p.3 line 9
			sss->b   = b;
			sss->n_b = n_b/2;
			SSS_query(sss);
		}
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Performs 'm'ost 's'ignificant 'b'it "less than" comparison.
 *  \return `true` if x < y; `false` otherwise.
 *
 *  [Chan] (p.4) proved: msb(x) < msb(y) iff ( x < y && x < (x^y) ).
 *
 *  <!-- References: -->
 *  [Chan]: ann/Chan2006_A_Minimalist's_Implementation_of_an_ANN_Algorithm_in_Fixed_Dimensions.pdf
 */
inline static bool less_msb
	(Index x, ///< 1st input.
	 Index y  ///< 2nd input.
	);

static void print_nodes (const ptrdiff_t n, const struct Node*const nodes, const char*const name)
{
	printf("Nodes (%s):\n",name);
	for (int i = 0; i < n; ++i) {
		const Index* xyz_i = nodes[i].xyz;
		for (int j = 0; j < DIM; ++j)
			printf(" %19d",xyz_i[j]);
		printf("\n");
	}
	printf("\n");
}

static const struct const_Matrix_d* constructor_xyz_min_max (const struct Input_ANN*const ann_i)
{
	struct Matrix_d*const xyz_mm = constructor_empty_Matrix_d('R',4,DIM); // returned
	double*const data_mm = xyz_mm->data;
	for (int d = 0; d < DIM; ++d) {
		data_mm[0*DIM+d] = DBL_MAX;
		data_mm[1*DIM+d] = DBL_MIN;
	}

	for (int i = 0; i < 2; ++i) {
		const struct const_Matrix_d*const b = ( i == 0 ? ann_i->nodes_b : ann_i->nodes_s );
		const ptrdiff_t n_b = b->ext_0;
		for (int n = 0; n < n_b; ++n) {
			const double*const data_n = get_row_const_Matrix_d(n,b);
			for (int d = 0; d < DIM; ++d) {
				if (data_n[d] < data_mm[0*DIM+d])
					data_mm[0*DIM+d] = data_n[d];
				if (data_n[d] > data_mm[1*DIM+d])
					data_mm[1*DIM+d] = data_n[d];
			}
		}
	}

	for (int d = 0; d < DIM; ++d) {
		data_mm[2*DIM+d] = 0.5*(data_mm[1*DIM+d] + data_mm[0*DIM+d]);
		data_mm[3*DIM+d] =      data_mm[1*DIM+d] - data_mm[0*DIM+d];
	}
	return (struct const_Matrix_d*) xyz_mm;
}

static int cmp_shuffle (Index* p, Index* q)
{
	int j = 0;
	Index x = 0;
	for (int k = 0; k < DIM; k++) {
		const Index y = (p[k]+shift)^(q[k]+shift);
		if (less_msb(x,y)) {
			j = k;
			x = y;
		}
	}
	return p[j]-q[j];
}

static void compute_distance_and_update (const ptrdiff_t n, struct SSS_c*const sss)
{
	struct Node*const p = &(sss->b[n]),
	           *const q = sss->s;

	Real z = 0;
	for (int j = 0; j < DIM; j++)
		z += POW2_R(p->xyz[j]-q->xyz[j]);

	Real*const r2 = &sss->r2;
	if (z < *r2) {
		sss->ind_ann = p->index;
		sss->ann     = p;
		*r2 = z;
		const Real r = sqrt(z); // Chan p.3 line 2
		for (int j = 0; j < DIM; j++) {
			// l == q^{s-[r]} (p.3, line 9)
			sss->l.xyz[j] = ( (q->xyz[j]-r > INDEX_MIN) ? (q->xyz[j]-(Index)ceil(r)) : INDEX_MIN );

			// u == q^{s+[r]} (p.3, line 6)
			sss->u.xyz[j] = ( (q->xyz[j]+r < INDEX_MAX) ? (q->xyz[j]+(Index)ceil(r)) : INDEX_MAX );
		}
	}
}

Real compute_r2_to_box (const struct SSS_c*const sss)
{
	const ptrdiff_t n = sss->n_b;
	const Index*const a_xyz = sss->b[0].xyz,
	           *const b_xyz = sss->b[n-1].xyz,
	           *const s_xyz = sss->s->xyz;

	int x = 0;
	for (int j = 0; j < DIM; ++j) {
		const Index y = (a_xyz[j]+shift)^(b_xyz[j]+shift);
		if (less_msb(x,y))
			x = y;
	}
	int i = -1;
	frexp(x,&i); // Extract most significant bit

	Real z = 0;
	for (int j = 0; j < DIM; ++j) {
		const Index x = ( (a_xyz[j]+shift) >>i ) << i,
		            y = x + (1 << i);
		if      (s_xyz[j]+shift < x)
			z += POW2_R(s_xyz[j]+shift-x);
		else if (s_xyz[j]+shift > y)
			z += POW2_R(s_xyz[j]+shift-y);
	}
	return z;
}

// Level 2 ********************************************************************************************************** //

inline static bool less_msb (Index x, Index y)
{
	return x < y && x < (x^y);
}
