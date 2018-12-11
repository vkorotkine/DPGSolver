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

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"

#include "matrix.h"
#include "vector.h"

#include "approximate_nearest_neighbor.h"

// Static function declarations ************************************************************************************* //

/// The shift to use. Global in this file as it must be used by \ref cmp_shuffle passed to qsort.
static Index shift = -INDEX_MAX;

/// \brief Set the value of the shift used for the ANN algorithm.
static void set_shift ( );

#define BASE_SEED 31415926535 ///< Random number base seed.

#define POW2_R(x) (((Real) (x))*((Real) (x))) ///< 'R'eal pow(x,2).

/// \brief Container for the converging node list and current search node.
struct SSS_c {
	ptrdiff_t n_b;      ///< The number of background nodes.
	struct Node_ANN* b; ///< The background nodes.
	struct Node_ANN* s; ///< The pointer to the search node.

	struct Node_ANN l; ///< Lower bounding node.
	struct Node_ANN u; ///< Upper bounding node.

	// ANN-related parameters
	Real r2; ///< The minimum squared Euclidian distance from the search node to computed background nodes.
	struct Node_ANN* ann; ///< The pointer to approximate nearest neighbor node.
};

/** \brief Constructor for a \ref SSS_ANN container.
 *  \return See brief. */
static const struct SSS_ANN* constructor_SSS_ANN
	(const struct Input_ANN*const ann_i ///< Standard.
	);

/// \brief Destructor for a \ref SSS_ANN container.
static void destructor_SSS_ANN
	(const struct SSS_ANN*const sss ///< Standard.
	);

/** \brief Perform the binary search for the input node in the list of background nodes.
 *
 *  ANN-related parameters of \ref SSS_c are set based on the current best guess before returning.
 */
static void SSS_query
	(struct SSS_c*const sss ///< Standard.
	);

/// \brief Print an array of \ref Node_ANN\*s.
static void print_nodes
	(const ptrdiff_t n,                 ///< Array length.
	 const struct Node_ANN*const nodes, ///< The array.
	 const char*const name              ///< Name to display.
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

/// \brief Copy the source into the destination \ref Node_ANN\*.
static void copy_Node_ANN
	(const struct Node_ANN*const src, ///< The source node.
	 struct Node_ANN*const dest       ///< The destination node.
	);

/// \brief Delete the node having index `ind_r` from the input array.
static void delete_Node_ANN
	(const ptrdiff_t ext_0,     ///< The length of the source array of nodes.
	 struct Node_ANN*const src, ///< The source array of nodes.
	 const ptrdiff_t ind_d      ///< The index of the node to delete.
	);

/** \brief Return the index in the input \ref Node_ANN array which has the input index.
 *  \return See brief. */
ptrdiff_t get_ind_Node_ANN_index
	(const ptrdiff_t index,          ///< The value of the index for which to search.
	 const ptrdiff_t ext_0,          ///< The length of the ANN node array.
	 const struct Node_ANN*const src ///< The ANN node array.
	);

// Interface functions ********************************************************************************************** //

const struct const_Vector_i* constructor_ann_indices (const struct Input_ANN*const ann_i)
{
	set_shift();
	const struct SSS_ANN*const sss = constructor_SSS_ANN(ann_i); // destructed
	const struct const_Vector_i*const ann_indices = constructor_ann_indices_from_sss(sss); // returned
	destructor_SSS_ANN(sss);

	return ann_indices;
}

const struct const_Vector_i* constructor_ann_indices_from_sss (const struct SSS_ANN*const sss)
{
	set_shift();
	const ptrdiff_t n_b = sss->ext_b,
	                n_s = sss->ext_s;

	struct Vector_i* ann_indices = constructor_empty_Vector_i(n_s); // returned
	for (int n = 0; n < n_s; ++n) {
		struct SSS_c sss_c =
			{ .n_b = n_b,
			  .b   = (struct Node_ANN*) sss->b,
			  .s   = (struct Node_ANN*) &(sss->s[n]),
			  .r2  = REAL_MAX,
			};
		SSS_query(&sss_c);
		ann_indices->data[n] = sss_c.ann->index;
		/* if (n == 0) { */
		/* 	printf("\n\nIndex (should be zero): %d\n\n\n",sss_c.ann->index); */
		/* 	EXIT; */
		/* } */
	}
	return (struct const_Vector_i*) ann_indices;
}

void destructor_Input_ANN (struct Input_ANN*const ann_info)
{
	destructor_const_Matrix_d(ann_info->nodes_b);
	destructor_const_Matrix_d(ann_info->nodes_s);
	free(ann_info);
}

void constructor_SSS_ANN_xyz (const struct Input_ANN*const ann_i, struct SSS_ANN*const sss)
{
	// Note: The same scaling must be used for all directions for the algorithm to work.
	struct Matrix_d*const xyz_mm = constructor_empty_Matrix_d('R',4,1); // keep
	double*const data_mm = xyz_mm->data;
	data_mm[0] = DBL_MAX;
	data_mm[1] = DBL_MIN;

	for (int i = 0; i < 2; ++i) {
		const struct const_Matrix_d*const b = ( i == 0 ? ann_i->nodes_b : ann_i->nodes_s );
		if (b == NULL)
			continue;
		assert(b->layout == 'R');
		const ptrdiff_t n_b = b->ext_0;
		for (int n = 0; n < n_b; ++n) {
			const double*const data_n = get_row_const_Matrix_d(n,b);
			for (int d = 0; d < DIM; ++d) {
				if (data_n[d] < data_mm[0])
					data_mm[0] = data_n[d];
				if (data_n[d] > data_mm[1])
					data_mm[1] = data_n[d];
			}
		}
	}
	data_mm[2] = 0.5*(data_mm[1] + data_mm[0]);
	data_mm[3] =      data_mm[1] - data_mm[0];

	sss->xyz_min_max = (struct const_Matrix_d*) xyz_mm;
}

void destructor_SSS_ANN_xyz (const struct SSS_ANN*const sss)
{
	destructor_const_Matrix_d(sss->xyz_min_max);
}

void constructor_SSS_ANN_b (const struct Input_ANN*const ann_i, struct SSS_ANN*const sss)
{
	enum { display = 0, };

	const ptrdiff_t n_b = ann_i->nodes_b->ext_0;

	sss->ext_b = n_b;
	struct Node_ANN*const b = calloc((size_t)n_b,sizeof *(sss->b)); // moved

	assert(ann_i->nodes_b != NULL);
	const struct const_Matrix_d*const xyz_mm = sss->xyz_min_max;
	const double*const xyz_avg = get_row_const_Matrix_d(2,xyz_mm),
	            *const xyz_h   = get_row_const_Matrix_d(3,xyz_mm);

	assert(ann_i->nodes_b->layout == 'R');
	const double* data_b = ann_i->nodes_b->data;
	for (int n = 0; n < n_b; ++n) {
		b[n].index = n;
		for (int d = 0; d < DIM; ++d) {
			const double rst = 2.0/xyz_h[0]*(*data_b-xyz_avg[0]);
			b[n].xyz[d] = (Index) (0.5*(((1.0-rst)*R_INDEX_MIN)+((1.0+rst)*R_INDEX_MAX)));
			++data_b;
		}
	}

	if (display)
		print_nodes((int)n_b,b,"background");
	sort_nodes_ANN(n_b,b);
	if (display)
		print_nodes((int)n_b,b,"background - sorted");

	sss->b = b; // keep
}

void destructor_SSS_ANN_b (const struct SSS_ANN*const sss)
{
	free((void*)sss->b);
}

void constructor_SSS_ANN_s (const struct Input_ANN*const ann_i, struct SSS_ANN*const sss)
{
	enum { display = 0, };

	const ptrdiff_t n_s = ann_i->nodes_s->ext_0;

	sss->ext_s = n_s;
	struct Node_ANN*const s = calloc((size_t)n_s,sizeof *(sss->s)); // moved

	assert(ann_i->nodes_s != NULL);
	const struct const_Matrix_d*const xyz_mm = sss->xyz_min_max;
	const double*const xyz_avg = get_row_const_Matrix_d(2,xyz_mm),
	            *const xyz_h   = get_row_const_Matrix_d(3,xyz_mm);

	assert(ann_i->nodes_s->layout == 'R');
	const double* data_s = ann_i->nodes_s->data;
	for (int n = 0; n < n_s; ++n) {
		s[n].index = n;
		for (int d = 0; d < DIM; ++d) {
			const double rst = 2.0/xyz_h[0]*(*data_s-xyz_avg[0]);
			s[n].xyz[d] = (Index) (0.5*((1.0-rst)*R_INDEX_MIN+(1.0+rst)*R_INDEX_MAX));
			++data_s;
		}
	}

	if (display)
		print_nodes((int)n_s,s,"search");

	sss->s = s; // keep
}

void destructor_SSS_ANN_s (const struct SSS_ANN*const sss)
{
	free((void*)sss->s);
}

void sort_nodes_ANN (const ptrdiff_t n_n, struct Node_ANN*const nodes)
{
	qsort((void*)nodes,(size_t)n_n,sizeof *nodes,(int (*)(const void *, const void *))cmp_shuffle);
}

const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_ANN (const struct const_Matrix_d*const nodes_i)
{
	set_shift();

	const ptrdiff_t n_n = nodes_i->ext_0;
	struct Vector_i*const inds_sorted = constructor_empty_Vector_i(n_n); // moved
	inds_sorted->data[0] = 0;

	struct Matrix_d*const nodes_b = constructor_empty_Matrix_d(nodes_i->layout,nodes_i->ext_0,nodes_i->ext_1); // moved
	struct Vector_d nodes_b_V = { .ext_0 = DIM, .data = NULL, };
	nodes_b_V.data = get_row_Matrix_d(0,nodes_b);
	set_to_data_Vector_d(&nodes_b_V,get_row_const_Matrix_d(0,nodes_i));


	struct SSS_ANN*const sss = calloc(1,sizeof *sss); // free

	struct Input_ANN ann_i = { .nodes_b = nodes_i, .nodes_s = NULL, };

	constructor_SSS_ANN_xyz(&ann_i,sss); // destructed
	constructor_SSS_ANN_b(&ann_i,sss);   // destructed
	const struct Node_ANN*const b0 = sss->b;
	struct Node_ANN s;
	sss->ext_s = 1;

	struct Node_ANN*const b_mutable = (struct Node_ANN*) sss->b;
	const ptrdiff_t ind_0 = get_ind_Node_ANN_index(0,sss->ext_b,sss->b);
	copy_Node_ANN(&sss->b[ind_0],&s);
	delete_Node_ANN(sss->ext_b,b_mutable,ind_0);
	sss->ext_b -= 1;

	for (int ind_s = 1; ind_s < n_n; ++ind_s) {
		struct SSS_c sss_c =
			{ .n_b = sss->ext_b,
			  .b   = (struct Node_ANN*) sss->b,
			  .s   = (struct Node_ANN*) &s,
			  .r2  = REAL_MAX,
			};
		SSS_query(&sss_c);
		const int ind_n = sss_c.ann->index;
		const ptrdiff_t ind_b = sss_c.ann-sss->b;
		inds_sorted->data[ind_s] = ind_n;

		copy_Node_ANN(&sss->b[ind_b],&s);
		delete_Node_ANN(sss->ext_b,b_mutable,ind_b);
		sss->ext_b -= 1;

		nodes_b_V.data = get_row_Matrix_d(ind_s,nodes_b);
		set_to_data_Vector_d(&nodes_b_V,get_row_const_Matrix_d(ind_n,nodes_i));
	}

	sss->ext_b = n_n;
	sss->b     = b0;

	destructor_SSS_ANN_xyz(sss);
	destructor_SSS_ANN_b(sss);
	free(sss);

	struct Nodes_Sorted_ANN*const nsa = calloc(1,sizeof *nsa); // returned
	nsa->nodes   = (struct const_Matrix_d*) nodes_b;     // keep
	nsa->indices = (struct const_Vector_i*) inds_sorted; // keep

	return nsa;
}

const struct Nodes_Sorted_ANN* constructor_Nodes_Sorted_ANN_with_trans (struct Matrix_d*const nodes_i)
{
	const bool requires_transpose = ( nodes_i->layout == 'R' ? false : true );
	if (requires_transpose)
		transpose_Matrix_d(nodes_i,true);

	const struct Nodes_Sorted_ANN*const ns = constructor_Nodes_Sorted_ANN((struct const_Matrix_d*)nodes_i);

	if (requires_transpose)
		transpose_Matrix_d(nodes_i,true);
	return ns;
}

void destructor_Nodes_Sorted_ANN (const struct Nodes_Sorted_ANN*const nsa)
{
	destructor_const_Matrix_d(nsa->nodes);
	destructor_const_Vector_i(nsa->indices);
	free((void*)nsa);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

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

static void set_shift ( )
{
	static bool needs_set = true;
	if (needs_set) {
		needs_set = false;
		srand48(BASE_SEED+DIM); // seed the random number generator.
		shift = (Index) (drand48()*R_INDEX_MAX);
	}
}

static const struct SSS_ANN* constructor_SSS_ANN (const struct Input_ANN*const ann_i)
{
	struct SSS_ANN*const sss = calloc(1,sizeof *sss); // free

	constructor_SSS_ANN_xyz(ann_i,sss); // destructed
	constructor_SSS_ANN_b(ann_i,sss);   // destructed
	constructor_SSS_ANN_s(ann_i,sss);   // destructed

	return sss;
}

static void destructor_SSS_ANN (const struct SSS_ANN*const sss)
{
	destructor_SSS_ANN_xyz(sss);
	destructor_SSS_ANN_b(sss);
	destructor_SSS_ANN_s(sss);
	free((void*)sss);
}

static void SSS_query (struct SSS_c*const sss)
{
	const ptrdiff_t n_b = sss->n_b;
	if (n_b == 0)
		return;

	compute_distance_and_update(n_b/2,sss);
	const Real r2 = sss->r2;

	if (n_b == 1 || compute_r2_to_box(sss)*POW2_R(1+EPS_ANN) > r2)
		return;

	struct Node_ANN*const b = sss->b,
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

static void print_nodes (const ptrdiff_t n, const struct Node_ANN*const nodes, const char*const name)
{
	printf("Nodes (%s):\n",name);
	for (int i = 0; i < n; ++i) {
		const Index* xyz_i = nodes[i].xyz;
		printf("%4d",nodes[i].index);
		for (int j = 0; j < DIM; ++j)
#if USE_SINGLE == true
			printf(" %19d",xyz_i[j]);
#else
			printf(" %19lld",xyz_i[j]);
#endif
		printf("\n");
	}
	printf("\n");
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
	if (p[j]-q[j] < 0)
		return -1;
	else if (p[j]-q[j] > 0)
		return 1;
	return 0;
}

static void copy_Node_ANN (const struct Node_ANN*const src, struct Node_ANN*const dest)
{
	dest->index = src->index;
	for (int d = 0; d < DIM; ++d)
		dest->xyz[d] = src->xyz[d];
}

static void delete_Node_ANN (const ptrdiff_t ext_0, struct Node_ANN*const src, const ptrdiff_t ind_d)
{
	memmove(&src[ind_d],&src[ind_d+1],(size_t)(ext_0-ind_d-1) * sizeof *src);
}

ptrdiff_t get_ind_Node_ANN_index (const ptrdiff_t index, const ptrdiff_t ext_0, const struct Node_ANN*const src)
{
	for (int i = 0; i < ext_0; ++i) {
		if (src[i].index == index)
			return i;
	}
	EXIT_ERROR("Did not find the index: %td\n",index);
}

// Level 1 ********************************************************************************************************** //

static void compute_distance_and_update (const ptrdiff_t n, struct SSS_c*const sss)
{
	struct Node_ANN*const p = &(sss->b[n]),
	               *const q = sss->s;

	Real z = 0;
	for (int j = 0; j < DIM; j++)
		z += POW2_R(p->xyz[j]-q->xyz[j]);

	Real*const r2 = &sss->r2;
	if (z < *r2) {
		sss->ann = p;
		*r2 = z;
		const Real r = sqrt(z); // Chan p.3 line 2
		for (int j = 0; j < DIM; j++) {
			// l == q^{s-[r]} (p.3, line 9)
			sss->l.xyz[j] = ( ((Real)q->xyz[j]-r > R_INDEX_MIN) ? (q->xyz[j]-(Index)ceil(r)) : INDEX_MIN );

			// u == q^{s+[r]} (p.3, line 6)
			sss->u.xyz[j] = ( ((Real)q->xyz[j]+r < R_INDEX_MAX) ? (q->xyz[j]+(Index)ceil(r)) : INDEX_MAX );
		}
	}
}

Real compute_r2_to_box (const struct SSS_c*const sss)
{
	const ptrdiff_t n = sss->n_b;
	const Index*const a_xyz = sss->b[0].xyz,
	           *const b_xyz = sss->b[n-1].xyz,
	           *const s_xyz = sss->s->xyz;

	Index x = 0;
	for (int j = 0; j < DIM; ++j) {
		const Index y = (a_xyz[j]+shift)^(b_xyz[j]+shift);
		if (less_msb(x,y))
			x = y;
	}
	int i = -1;
	frexp((Real)x,&i); // Extract most significant bit

	Real z = 0;
	for (int j = 0; j < DIM; ++j) {
		const Index x = ( (a_xyz[j]+shift) >> i ) << i,
		            y = x + ((Index)1 << i);
		if      (s_xyz[j]+shift < x)
			z += POW2_R(s_xyz[j]+shift-x);
		else if (s_xyz[j]+shift > y)
			z += POW2_R(s_xyz[j]+shift-y);
	}
	return z;
}

inline static bool less_msb (Index x, Index y)
{
	return x < y && x < (x^y);
}
