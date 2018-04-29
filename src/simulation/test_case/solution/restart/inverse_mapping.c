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

#include "inverse_mapping.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_elements.h"
#include "definitions_tol.h"

#include "matrix.h"

#include "math_functions.h"
#include "nodes.h"

// Static function declarations ************************************************************************************* //

struct Newton;

/** \brief Pointer \ref Newton data computing functions.
 *  \param ve  The matrix of vertices.
 *  \param xyz The xyz data.
 *  \param rst The rst coordinates (standard reference element).
 */
typedef struct Newton (*get_Newton_fptr)
	(const struct const_Matrix_d*const ve,
	 const double*const xyz,
	 const double*const rst
	);

/** \brief Pointer \ref Newton update term computing functions.
 *  \param newton \ref Newton.
 */
typedef const double* (*get_Newton_term_fptr)
	(const struct Newton*const newton
	);

/// \brief Container for data related to Newton's method.
struct Newton {
	double res[DIM];     ///< The residual terms.
	double jac[DIM*DIM]; ///< The residual jacobian terms.
};

/// \brief Container for functionality related to the computation of the inverse mapping.
struct Inv_Map {
	int dim;    ///< The dimension of the nodes.
	int type;   ///< \ref Element::type.
	int s_type; ///< \ref Element::s_type.

	double guess; ///< Value for the initial guess for the Newton algorithm (close to the centroid is best).

	get_Newton_fptr get_Newton;           ///< Standard.
	get_Newton_term_fptr get_Newton_term; ///< Standard.
};

/** \brief Constructor for a \ref Inv_Map container.
 *  \return See brief. */
static const struct Inv_Map* constructor_Inv_Map
	(const int e_type ///< \ref Element::type.
	);

/// \brief Destructor for a \ref Inv_Map container.
static void destructor_Inv_Map
	(const struct Inv_Map*const im ///< Standard.
	);

/** \brief Version of \ref get_Newton_fptr for LINEs.
 *  \return See brief. */
static struct Newton get_Newton_line
	(const struct const_Matrix_d*const ve, ///< See brief.
	 const double*const xyz,               ///< See brief.
	 const double*const rst                ///< See brief.
	);

/** \brief Version of \ref get_Newton_fptr for TRIs.
 *  \return See brief. */
static struct Newton get_Newton_tri
	(const struct const_Matrix_d*const ve, ///< See brief.
	 const double*const xyz,               ///< See brief.
	 const double*const rst                ///< See brief.
	);

/** \brief Version of \ref get_Newton_fptr for QUADs.
 *  \return See brief. */
static struct Newton get_Newton_quad
	(const struct const_Matrix_d*const ve, ///< See brief.
	 const double*const xyz,               ///< See brief.
	 const double*const rst                ///< See brief.
	);

/** \brief Version of \ref get_Newton_term_fptr for 1 dimensional term.
 *  \return See brief. */
static const double* get_Newton_term_1d
	(const struct Newton*const newton ///< See brief.
	);

/** \brief Version of \ref get_Newton_term_fptr for 2 dimensional term.
 *  \return See brief. */
static const double* get_Newton_term_2d
	(const struct Newton*const newton ///< See brief.
	);

/** \brief Constructor for reference rst coordinates from standard reference element reference coordinates.
 *  \return See brief.
 *
 *  The reference element used in this code can be identified according to the vertices in
 *  \ref constructor_const_Nodes_vertices. The standard reference element has vertices between in [0,1], generally
 *  having vertices determined from the encompasing HEX.
 */
static const struct const_Matrix_d* constructor_rst_ref_from_rst_std
	(const struct Inv_Map*const im,            ///< Standard.
	 const struct const_Matrix_d*const rst_std ///< The standard reference element coordinates.
	);

// Interface functions ********************************************************************************************** //

struct Matrix_d* constructor_inverse_mapping_mutable
	(const int e_type, const struct const_Matrix_d*const xyz_ve, const struct const_Matrix_d*const xyz)
{
	/** Newton's method is used here to converge to the rst coordinates which satisfy: xyz = xyz_ve*phi(rst).
	 *
	 *  In the case of affine elements (simplices), the algorithm converges in a single iteration as the function
	 *  above is linear in rst. In other cases, the convergence is quadratic and it is not necessary to treat
	 *  special cases (such as for parallelogram vertices).
	 *
	 *  Please consult [this][SO_inv_map] stack overflow answer for additional details.
	 *
	 *  <!-- References: -->
	 *  [SO_inv_map]: https://stackoverflow.com/a/18332009/5983549
	 */
	enum { count_max = 10, }; // The maximum number of Newton steps before assuming no convergence.

	assert(xyz_ve->layout == 'R');
	assert(xyz->layout == 'R');

	const ptrdiff_t n_n = xyz->ext_0,
	                dim = xyz->ext_1;

	const struct Inv_Map*const im = constructor_Inv_Map(e_type); // destructed
	assert(dim == im->dim);

	struct Matrix_d*const rst = constructor_empty_Matrix_d('R',n_n,dim); // destructed
	set_to_value_Matrix_d(rst,im->guess);

	for (int n = 0; n < n_n; ++n) {
		const double*const data_xyz = get_row_const_Matrix_d(n,xyz);
		double*const data_rst       = get_row_Matrix_d(n,rst);

		int count = 0;
		for ( ; count < count_max; ++count) {
			struct Newton newton = im->get_Newton(xyz_ve,data_xyz,data_rst);
			if (maximum_abs_d(newton.res,dim) < EPS)
				break;

			const double*const drst = im->get_Newton_term(&newton);
			for (int d = 0; d < dim; ++d)
				data_rst[d] -= drst[d];
		}
		assert(count < count_max); // Not yet converged.
	}

	transpose_Matrix_d(rst,true);
	struct Matrix_d*const rst_ref =
		(struct Matrix_d*) constructor_rst_ref_from_rst_std(im,(struct const_Matrix_d*)rst); // returned
	destructor_Matrix_d(rst);
	destructor_Inv_Map(im);

	return rst_ref;
}

const struct const_Matrix_d* constructor_inverse_mapping
	(const int e_type, const struct const_Matrix_d*const xyz_ve, const struct const_Matrix_d*const xyz)
{
	return (struct const_Matrix_d*) constructor_inverse_mapping_mutable(e_type,xyz_ve,xyz);
}

const struct const_Matrix_d* constructor_basis_std_p1 (const int e_type, const struct const_Matrix_d*const rst_std)
{
	const ptrdiff_t n_n = rst_std->ext_0,
	                dim = rst_std->ext_1;

	const double*const r = get_col_const_Matrix_d(0,rst_std),
	            *const s = ( dim > 1 ? get_col_const_Matrix_d(1,rst_std) : NULL),
	            *const t = ( dim > 2 ? get_col_const_Matrix_d(2,rst_std) : NULL);

	struct Matrix_d* phi_std = NULL;
	switch (e_type) {
	case LINE: {
		assert(dim == 1);
		const ptrdiff_t n_b = dim+1;
		phi_std = constructor_empty_Matrix_d('R',n_n,n_b); // returned
		for (int n = 0; n < n_n; ++n) {
			double* data_phi = get_row_Matrix_d(n,phi_std);
			*data_phi++ = 1.0-r[n];
			*data_phi++ = r[n];
		}
		break;
	} case TRI: {
		assert(dim == 2);
		const ptrdiff_t n_b = dim+1;
		phi_std = constructor_empty_Matrix_d('R',n_n,n_b); // returned
		for (int n = 0; n < n_n; ++n) {
			double* data_phi = get_row_Matrix_d(n,phi_std);
			*data_phi++ = 1.0-r[n]-s[n];
			*data_phi++ = r[n];
			*data_phi++ = s[n];
		}
		break;
	} case QUAD: {
		assert(dim == 2);
		const ptrdiff_t n_b = 4;
		phi_std = constructor_empty_Matrix_d('R',n_n,n_b); // returned
		for (int n = 0; n < n_n; ++n) {
			double* data_phi = get_row_Matrix_d(n,phi_std);
			*data_phi++ = (1.0-r[n])*(1.0-s[n]);
			*data_phi++ = (    r[n])*(1.0-s[n]);
			*data_phi++ = (1.0-r[n])*(    s[n]);
			*data_phi++ = (    r[n])*(    s[n]);
		}
		break;
	} case TET: {
		assert(dim == 1);
		const ptrdiff_t n_b = dim+1;
		phi_std = constructor_empty_Matrix_d('R',n_n,n_b); // returned
EXIT_ADD_SUPPORT; UNUSED(t);
	} default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
	return (struct const_Matrix_d*) phi_std;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct Inv_Map* constructor_Inv_Map (const int e_type)
{
	struct Inv_Map*const im = calloc(1,sizeof *im); // free

	im->type = e_type;
	switch (e_type) {
	case LINE:
		im->dim    = 1;
		im->s_type = ST_TP;
		im->guess  = 1.0/2.0;
		im->get_Newton      = get_Newton_line;
		im->get_Newton_term = get_Newton_term_1d;
		break;
	case TRI:
		im->dim    = 2;
		im->s_type = ST_SI;
		im->guess  = 1.0/3.0;
		im->get_Newton      = get_Newton_tri;
		im->get_Newton_term = get_Newton_term_2d;
		break;
	case QUAD:
		im->dim    = 2;
		im->s_type = ST_TP;
		im->guess  = 1.0/2.0;
		im->get_Newton      = get_Newton_quad;
		im->get_Newton_term = get_Newton_term_2d;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",e_type);
		break;
	}
	return im;
}

static void destructor_Inv_Map (const struct Inv_Map*const im)
{
	free((void*)im);
}

static struct Newton get_Newton_line
	(const struct const_Matrix_d*const ve, const double*const xyz, const double*const rst)
{
	assert(ve->layout == 'R');

	enum { dim = 1, };
	assert(ve->ext_1 == dim);

	const double*const p  = xyz,
	            *const p0 = ve->data,
	            *const p1 = p0+dim;
	const double r = rst[0];

	static struct Newton newton;
	double*const res = newton.res,
	      *const jac = newton.jac;

	int ind = 0;
	for (int i = 0; i < dim; ++i) {
		res[ind++] = p[i] - ( p0[i]*(1-r) + p1[i]*(r) );
	}
	ind = 0;
	for (int i = 0; i < dim; ++i) {
//		jac[ind++] = - ( p0[i]*(0-1) + p1[i]*(1) );
		jac[ind++] = p0[i] - p1[i];
	}
	return newton;
}

static struct Newton get_Newton_tri
	(const struct const_Matrix_d*const ve, const double*const xyz, const double*const rst)
{
	assert(ve->layout == 'R');

	enum { dim = 2, };
	assert(ve->ext_1 == dim);

	const double*const p  = xyz,
	            *const p0 = ve->data,
	            *const p1 = p0+dim,
	            *const p2 = p1+dim;
	const double r = rst[0],
	             s = rst[1];

	static struct Newton newton;
	double*const res = newton.res,
	      *const jac = newton.jac;

	int ind = 0;
	for (int i = 0; i < dim; ++i) {
		res[ind++] = p[i] - ( p0[i]*(1-r-s) + p1[i]*(r) + p2[i]*(s) );
	}
	ind = 0;
	for (int i = 0; i < dim; ++i) {
//		jac[ind++] = - ( p0[i]*(0-1-0) + p1[i]*(1) + p2[i]*(0) );
//		jac[ind++] = - ( p0[i]*(0-0-1) + p1[i]*(0) + p2[i]*(1) );
		jac[ind++] = p0[i] - p1[i];
		jac[ind++] = p0[i] - p2[i];
	}
	return newton;
}

static struct Newton get_Newton_quad
	(const struct const_Matrix_d*const ve, const double*const xyz, const double*const rst)
{
	assert(ve->layout == 'R');

	enum { dim = 2, };
	assert(ve->ext_1 == dim);

	const double*const p  = xyz,
	            *const p0 = ve->data,
	            *const p1 = p0+dim,
	            *const p2 = p1+dim,
	            *const p3 = p2+dim;
	const double r = rst[0],
	             s = rst[1];

	static struct Newton newton;
	double*const res = newton.res,
	      *const jac = newton.jac;

	int ind = 0;
	for (int i = 0; i < dim; ++i) {
		res[ind++] = p[i] - ( p0[i]*(1-r)*(1-s) + p1[i]*(r)*(1-s) + p2[i]*(1-r)*(s) + p3[i]*(r)*(s) );
	}
	ind = 0;
	for (int i = 0; i < dim; ++i) {
//		jac[ind++] = - ( p0[i]*(0-1)*(1-s) + p1[i]*(1)*(1-s) + p2[i]*(0-1)*(s) + p3[i]*(1)*(s) );
//		jac[ind++] = - ( p0[i]*(1-r)*(0-1) + p1[i]*(r)*(0-1) + p2[i]*(1-r)*(1) + p3[i]*(r)*(1) );
		jac[ind++] = p0[i]*(1-s) - p1[i]*(1-s) + p2[i]*(s) - p3[i]*(s);
		jac[ind++] = p0[i]*(1-r) + p1[i]*(r) - p2[i]*(1-r) - p3[i]*(r);
	}
	return newton;
}

static const double* get_Newton_term_1d (const struct Newton*const newton)
{
	enum { dim = 1, };
	const double*const res = newton->res,
	            *const jac = newton->jac;

	const double inv_det_jac = 1.0/(jac[0]);
	const double inv_jac[] = {  1.0, };

	static double drst[dim] = { 0.0, };
	for (int i = 0; i < dim; ++i) {
		drst[i] = 0.0;
		for (int j = 0; j < dim; ++j)
			drst[i] += inv_jac[i*dim+j]*res[j];
		drst[i] *= inv_det_jac;
	}
	return drst;
}

static const double* get_Newton_term_2d (const struct Newton*const newton)
{
	enum { dim = 2, };
	const double*const res = newton->res,
	            *const jac = newton->jac;

	const double inv_det_jac = 1.0/(jac[0]*jac[3]-jac[1]*jac[2]);
	const double inv_jac[] = {  jac[3], -jac[1], -jac[2], jac[0], };

	static double drst[dim] = { 0.0, };
	for (int i = 0; i < dim; ++i) {
		drst[i] = 0.0;
		for (int j = 0; j < dim; ++j)
			drst[i] += inv_jac[i*dim+j]*res[j];
		drst[i] *= inv_det_jac;
	}
	return drst;
}

static const struct const_Matrix_d* constructor_rst_ref_from_rst_std
	(const struct Inv_Map*const im, const struct const_Matrix_d*const rst_std)
{
	const struct const_Matrix_d*const phi_rst = constructor_basis_std_p1(im->type,rst_std);             // dest.
	const struct const_Nodes*const ve_nodes   = constructor_const_Nodes_vertices(im->dim,1,im->s_type); // dest.

	const struct const_Matrix_d*const rst_ref =
		constructor_mm_const_Matrix_d('N','N',1.0,phi_rst,ve_nodes->rst,'R'); // returned
	destructor_const_Matrix_d(phi_rst);
	destructor_const_Nodes(ve_nodes);

	return rst_ref;
}
