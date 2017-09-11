// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__cubature_h__INCLUDED
#define DPG__cubature_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions computing the reference element coordinate (and associated cubature
 *	       if relevant).
 */

#include <stdbool.h>

/// Container for reference element coordinates and cubature related information.
struct Cubature {
	int p,         ///< The order of the nodes.
	    node_type; ///< The index associated with the node type.

	struct Matrix_d* rst; ///< The reference element coordinates.

	bool has_weights;   ///< Flag for whether weights are included.
	struct Vector_d* w; ///< The cubature weights.
};

/** \brief Constructor for a \ref Cubature container of tensor-product type.
 *  \return Standard.
 *
 *  The implementation is based off of that of Shen et al. (section 3.3.2, \cite Shen2011), with implementation scripts
 *  provided on the [associated website][shen_web].
 *
 *  \todo Needs modernizing.
 *
 *  <!-- References: -->
 *  [shen_web]: http://www.ntu.edu.sg/home/lilian/book.htm
 */
struct Cubature* constructor_Cubature_TP
	(const int d,        ///< The dimension of the nodes.
	 const int p,        ///< Defined in \ref Cubature.
	 const int node_type ///< Defined in \ref Cubature.
	);

/*
typedef void (*cubature_tdef) (struct S_CUBATURE *const CUBDATA);

extern void cubature_TP  (struct S_CUBATURE *const CUBDATA);
extern void cubature_TRI (struct S_CUBATURE *const CUBDATA);
extern void cubature_TET (struct S_CUBATURE *const CUBDATA);
extern void cubature_PYR (struct S_CUBATURE *const CUBDATA);

extern void set_cubdata (struct S_CUBATURE *const CUBDATA, bool const return_w, bool const return_symm,
                         char const *const NodeType, unsigned int const d, unsigned int const P,
                         cubature_tdef cubature);
extern void set_from_cubdata (struct S_CUBATURE const *const CUBDATA, unsigned int *Nn, unsigned int *Ns, double **rst,
                              double **w, unsigned int **symms);

extern struct S_CUBATURE *cub_constructor (bool const return_w, bool const return_symm, char const *const NodeType,
                                           unsigned int const d, unsigned int const P, cubature_tdef cubature);
extern void cub_destructor (struct S_CUBATURE *CUBDATA);
*/

#endif // DPG__cubature_h__INCLUDED
