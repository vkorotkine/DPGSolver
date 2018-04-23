// Converting to c and cleaning up.

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <limits.h>
#include <math.h>

//#define SINGLE // single or double precision
#ifdef SINGLE
	#define Real      float
	#define Type      int
	#define Index     int
	#define INDEX_MAX INT_MAX
	#define REAL_MAX  FLT_MAX
#else
	#define Real      double
	#define Type      long long // Change to double after testing.
	#define Index     long long
	#define INDEX_MAX LLONG_MAX
	#define REAL_MAX  DBL_MAX
#endif

#define N_B 16
#define N_S 1
#define DIM 2


#define SQ(x) (((Real) (x))*((Real) (x)))

#define PRINT_FILELINE ({printf("\n\nPrinting at: FILE: %s, FUNCTION: %s (LINE: %d)\n\n",__FILE__,__func__,__LINE__);})
#define EXIT_MSG         ({ PRINT_FILELINE; fflush(stdout); abort(); })
#define EXIT_UNSUPPORTED ({printf("\n\nError: Unsupported.\n"), EXIT_MSG; })
#define EXIT_ADD_SUPPORT ({printf("\n\nError: Add support.\n"), EXIT_MSG; })
#define EXIT_ERROR(...)  ({printf("\n\nError: "); printf(__VA_ARGS__); printf("\n\n"); EXIT_MSG; })
#define UNUSED(x)       (void)(x)


// Static function declarations ************************************************************************************* //

Index shift;        ///< The shift to use. Global in this file as it must be used by cmp_shuffle passed to qsort.
#define EPS_ANN 0.0 ///< The accuracy tolerance for the approximate nearest neighbor search.

/// \brief Container for a node represented in binary.
struct Node {
	Index xyz[DIM]; ///< The coordinates.
	int index;      ///< The index.
};

/// \brief Container for 'S'hift-'S'huffle-'S'ort approximate nearest neighbor information.
struct SSS {
	struct Node p[N_B]; ///< The background nodes.
	struct Node q[N_S]; ///< The search nodes.
};

/// \brief Container for the converging node list and current search node.
struct SSS_c {
	int n_b;        ///< The number of background nodes.
	struct Node* p; ///< The background nodes.
	struct Node* q; ///< The pointer to the search node.

	struct Node a; ///< Lower bounding node.
	struct Node b; ///< Upper bounding node.

	Real r2;     ///< The minimum Euclidian distance from the search node to computed background nodes.
	int ind_ann; ///< The current index of the approximate nearest neighbor.
};

/** \brief Sort the nodes according to the shuffle order.
 *
 *  See Chan section 2.
 */
static void SSS_preprocess
	(struct SSS*const sss ///< \ref SSS.
	);

/** \brief Perform the binary search for the input node in the list of background nodes.
 *
 *  \ref SSS_c::ind_ann is set to the value of the index of the background node which is the approximate nearest
 *  neighbor after returning.
 */
static void SSS_query
	(struct SSS_c*const sss ///< \ref SSS_c.
	);

// Interface functions ********************************************************************************************** //

int main (int argc, char** argv)
{
	srand48(31415+N_B+N_S+DIM); // seed the random number generator.
	shift = (Index) (drand48()*INDEX_MAX);

	struct SSS sss;

	if (sizeof(Type) != sizeof(Index))
		EXIT_ERROR("Incompatible sizes: %zu != %zu.",sizeof(Type),sizeof(Index));

	const Type p_scale = 4;
	const Type p_data[N_B*DIM] = { 0,0, 0,1, 0,2, 0,3, 1,0, 1,1, 1,2, 1,3,
	                               2,0, 2,1, 2,2, 2,3, 3,0, 3,1, 3,2, 3,3, };
	for (int i = 0, ind = 0; i < N_B; i++) {
		sss.p[i].index  = i;
		for (int j = 0; j < DIM; j++) {
			sss.p[i].xyz[j] = (Index) p_scale*p_data[ind];
			++ind;
		}
	}

//	const Type q_data[N_S*DIM] = { 13.1, 8.9, };
	const Type q_data[N_S*DIM] = { 13, 9, };
	for (int i = 0, ind = 0; i < N_S; i++) {
	for (int j = 0; j < DIM; j++) {
		sss.q[i].xyz[j] = (Index) q_data[ind++];
	}}

	SSS_preprocess(&sss);

	for (int i = 0; i < N_S; i++) {
		struct SSS_c sss_c =
			{ .r2  = REAL_MAX,
			  .n_b = N_B,
		        .p   = sss.p,
		        .q   = &(sss.q[i]),
		      };
		SSS_query(&sss_c);
		const struct Node*const node_ann = &(sss_c.p[sss_c.ind_ann]);
		printf("%d %d\n",i,node_ann->index);
	}
EXIT_ADD_SUPPORT;

	return 0;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Shuffle order comparison function for nodes.
 *
 *  Chan (p.4) proved: p <~ q iff x_j <= y_j, j is the index of the 'm'ost 's'ignificant 'b'it of (x xor y) for j \in d.
 *  The <~ operator, called the shuffle order, is used to denote "less than" comparison of the shuffle of two points
 *  written in binary.
 */
static int cmp_shuffle
	(Index* p, ///< The 1st node (represented in binary).
	 Index* q  ///< The 2nd node (represented in binary).
	);

static void SSS_preprocess (struct SSS*const sss)
{
	qsort((void*)sss->p,N_B,sizeof(struct Node),(int (*)(const void *, const void *))cmp_shuffle);
}

// Level 1 ********************************************************************************************************** //

/** \brief Performs 'm'ost 's'ignificant 'b'it "less than" comparison.
 *
 *  Chan (p.4) proved: msb(x) < msb(y) iff ( x < y && x < (x^y) ).
 */
inline static bool less_msb (Index x, Index y) { return x < y && x < (x^y); }

/** \brief Compute the distance between the background node of index n and the search node and update
 *         \ref SSS_c::ind_ann if the current node is closer than all previous. */
static void compute_distance_and_update
	(const int n,           ///< The index of \ref SSS_c::p.
	 struct SSS_c*const sss ///< \ref SSS_c.
	);

/** \brief Return the squared distance from the search node to the current bounding box.
 *  \return See brief. */
static Real compute_r2_to_box
	(const struct SSS_c*const sss ///< \ref SSS_c.
	);

static int cmp_shuffle (Index* p, Index* q)
{
	int j = 0;
	Index x = 0,
	      y = 0;
	for (int k = 0; k < DIM; k++) {
		const Index y = (p[k]+shift)^(q[k]+shift);
		if (less_msb(x,y)) {
			j = k;
			x = y;
		}
	}
	return p[j]-q[j];
}

static void SSS_query (struct SSS_c*const sss)
{
	const int n_b = sss->n_b;
	if (n_b == 0)
		return;

	compute_distance_and_update(n_b/2,sss);
	const Real r2 = sss->r2;

//	printf("% 3d % .3e % .3e\n",n_b,dist_sq_to_box(q, P[0],P[n-1])*sq(1+eps),r_sq);
struct Node* p1 = sss->p;
	printf("% 3d %llu %llu %llu %llu % .3e % .3e\n",
	       n_b,p1[0].xyz[0],p1[0].xyz[1],p1[n_b-1].xyz[0],p1[n_b-1].xyz[1],compute_r2_to_box(sss)*SQ(1+EPS_ANN),r2);

	if (n_b == 1 || compute_r2_to_box(sss)*SQ(1+EPS_ANN))
		return;

	struct Node*const p = sss->p,
	           *const q = sss->q;

	if (cmp_shuffle(q->xyz,p[n_b/2].xyz) < 0) { // p.3 line 4
		sss->n_b = n_b/2; // Chan p.3 line 5 (binary search in lower half)
		SSS_query(sss);

		if (cmp_shuffle(sss->b.xyz,p[n_b/2].xyz) > 0) { // Chan p.3 line 6
			sss->p   = &p[n_b/2+1];
			sss->n_b = n_b-(n_b/2+1);
			SSS_query(sss);
		}
	} else {
		sss->p   = &p[n_b/2+1];   // p.3 line 8 (binary search in upper half)
		sss->n_b = n_b-(n_b/2+1);
		SSS_query(sss);

		if (cmp_shuffle(sss->a.xyz,p[n_b/2].xyz) < 0) { // Chan p.3 line 9
			sss->n_b = n_b/2;
			SSS_query(sss);
		}
	}

EXIT_ADD_SUPPORT; UNUSED(sss);
}

// Level 2 ********************************************************************************************************** //

static void compute_distance_and_update (const int n, struct SSS_c*const sss)
{
	struct Node*const p = &(sss->p[n]),
	           *const q = sss->q;

	Real z = 0;
	for (int j = 0; j < DIM; j++)
{
//printf("%llu %llu %llu\n",p->xyz[j],q->xyz[j],p->xyz[j]-q->xyz[j]);
		z += SQ(p->xyz[j]-q->xyz[j]);
}
//EXIT_UNSUPPORTED;

	Real*const r2 = &sss->r2;
//printf("% .3e % .3e\n",z,*r2);
	if (z < *r2) {
		sss->ind_ann = p->index;
		*r2 = z;
		const Real r = sqrt(z); // Chan p.3 line 2
		for (int j = 0; j < DIM; j++) {
			// a == q^{s-[r]} (p.3, line 9)
			sss->a.xyz[j] = (q->xyz[j]   > r)         ? (q->xyz[j]-(Index)ceil(r)) : 0;

			// b == q^{s+[r]} (p.3, line 6)
			sss->b.xyz[j] = (q->xyz[j]+r < INDEX_MAX) ? (q->xyz[j]+(Index)ceil(r)) : INDEX_MAX;
		}
	}
}

static Real compute_r2_to_box (const struct SSS_c*const sss)
{
	const int n = sss->n_b;
	const Index*const a_xyz = sss->p[0].xyz,
	           *const b_xyz = sss->p[n-1].xyz,
	           *const q_xyz = sss->q->xyz;

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
		if      (q_xyz[j]+shift < x)
			z += SQ(q_xyz[j]+shift-x);
		else if (q_xyz[j]+shift > y)
			z += SQ(q_xyz[j]+shift-y);
	}
	return z;
}
