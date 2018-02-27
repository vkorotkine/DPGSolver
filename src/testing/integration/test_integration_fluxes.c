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

#include <assert.h>
#include <string.h>
#include "petscsys.h"

#include "macros.h"
#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_test_integration.h"
#include "definitions_tol.h"

#include "test_base.h"
#include "test_integration.h"
#include "test_support.h"
#include "test_support_multiarray.h"

#include "multiarray.h"

#include "flux.h"
#include "simulation.h"
#include "test_case.h"

#include "complex_multiarray.h"

#include "test_complex_flux.h"
#include "test_complex_math_functions.h"
#include "test_complex_test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Constructor for the data members of the real \ref Flux_Input_T container.
static void constructor_Flux_Input_data_members
	(struct Flux_Input*const flux_i,   ///< \ref Flux_Input_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Destructor for the data members of the real \ref Flux_Input_T container.
static void destructor_Flux_Input_data_members
	(struct Flux_Input*const flux_i ///< \ref Flux_Input_T.
	);

/// \brief Constructor for the data members of the complex \ref Flux_Input_T container.
static void constructor_Flux_Input_c_data_members
	(struct Flux_Input_c*const flux_c_i, ///< \ref Flux_Input_T.
	 struct Flux_Input*const flux_i      ///< \ref Flux_Input_T.
	);

/// \brief Destructor for the data members of the complex \ref Flux_Input_T container.
static void destructor_Flux_Input_c_data_members
	(struct Flux_Input_c*const flux_c_i ///< \ref Flux_Input_T.
	);

/** \brief Constructor for the real \ref Flux_T computed using the complex step method.
 *  \return See brief. */
static struct Flux* constructor_Flux_cmplx_step
	(struct Flux_Input_c*const flux_i_c, ///< Container for the complex flux input data.
	 struct Flux_Input*const flux_i      ///< Container for the real flux input data.
	);

// Interface functions ********************************************************************************************** //

/** \test Performs integration testing for the fluxes (\ref test_integration_fluxes.c).
 *  \return 0 on success.
 *
 *  The flux linearizations are verified by comparing with the output when using the complex step method. Details of the
 *  complex step method can be found in Squire et al. \cite Squire1998 and Martins et al. \cite Martins2003.
 */
int main
	(int argc,   ///< Standard.
	 char** argv ///< Standard.
	)
{
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

	assert_condition_message(argc == 3,"Invalid number of input arguments");
	const char*const ctrl_name = argv[2];

	struct Integration_Test_Info*const int_test_info =
		constructor_Integration_Test_Info(ctrl_name); // destructed

	const int p          = int_test_info->p_ref[0],
	          ml         = int_test_info->ml[0],
	          adapt_type = int_test_info->adapt_type;
	destructor_Integration_Test_Info(int_test_info);
	assert(adapt_type == ADAPT_0);

	const char*const ctrl_name_curr = set_file_name_curr(adapt_type,p,ml,false,ctrl_name);

	const char type_rc[2] = { 'r', 'c', };
	struct Simulation* sim[2] = { NULL, NULL, };
	int ind_rc = -1;

	// real
	ind_rc = 0;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc]); // destructed
	((struct Test_Case*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'i';

	struct Flux_Input* flux_i = constructor_Flux_Input(sim[ind_rc]); // destructed
	constructor_Flux_Input_data_members(flux_i,sim[ind_rc]);         // destructed
	struct Flux* flux = constructor_Flux(flux_i);                    // destructed

	// complex
	ind_rc = 1;
	structor_simulation(&sim[ind_rc],'c',adapt_type,p,ml,0,0,ctrl_name_curr,type_rc[ind_rc]); // destructed
	convert_to_Test_Case_rc(sim[ind_rc],'c');
	((struct Test_Case_c*)sim[ind_rc]->test_case_rc->tc)->solver_method_curr = 'e';

	struct Flux_Input_c* flux_c_i = constructor_Flux_Input_c(sim[ind_rc]); // destructed
	constructor_Flux_Input_c_data_members(flux_c_i,flux_i);                // destructed
	struct Flux_c* flux_c = constructor_Flux_c(flux_c_i);                  // destructed
	destructor_Flux_c(flux_c);

	struct Test_Info test_info = { .n_warn = 0, };

	struct Flux* flux_cmplx_step = constructor_Flux_cmplx_step(flux_c_i,flux_i);
	convert_to_Test_Case_rc(sim[ind_rc],'r');

	const bool* compute_member = flux_i->compute_member;

	assert(!compute_member[4]);
	assert(!compute_member[5]);
	bool pass        = false;
	const double tol = EPS;
	const bool differences[] =
		{ compute_member[0] ? diff_const_Multiarray_d(flux->f,      flux_cmplx_step->f,      tol) : false,
		  compute_member[1] ? diff_const_Multiarray_d(flux->df_ds,  flux_cmplx_step->df_ds,  tol) : false,
		  compute_member[2] ? diff_const_Multiarray_d(flux->df_dg,  flux_cmplx_step->df_dg,  tol) : false,
		  compute_member[3] ? diff_const_Multiarray_d(flux->d2f_ds2,flux_cmplx_step->d2f_ds2,tol) : false,
		};
	const int len = COUNT_OF(differences);
	if (check_diff(len,differences,&pass)) {
		if (differences[0]) print_diff_const_Multiarray_d(flux->f,      flux_cmplx_step->f,      tol);
		if (differences[1]) print_diff_const_Multiarray_d(flux->df_ds,  flux_cmplx_step->df_ds,  tol);
		if (differences[2]) print_diff_const_Multiarray_d(flux->df_dg,  flux_cmplx_step->df_dg,  tol);
		if (differences[3]) print_diff_const_Multiarray_d(flux->d2f_ds2,flux_cmplx_step->d2f_ds2,tol);
	}
	expect_condition(pass,"flux_linearization");
	assert_condition(pass);

	destructor_Flux(flux_cmplx_step);

	// complex
	ind_rc = 1;
	destructor_Flux_Input_c_data_members(flux_c_i);
	destructor_Flux_Input_c(flux_c_i);
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc]);

	// real
	ind_rc = 0;
	destructor_Flux(flux);
	destructor_Flux_Input_data_members(flux_i);
	destructor_Flux_Input(flux_i);
	structor_simulation(&sim[ind_rc],'d',adapt_type,p,ml,0,0,NULL,type_rc[ind_rc]);

	output_warning_count(&test_info);

	PetscFinalize();
	OUTPUT_SUCCESS;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Constructor for the 's'olution data.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_s
	(const int index,                  ///< The index of the available solution variables.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Constructor for the solution 'g'radient data.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_g
	(const int index,                  ///< The index of the available solution variables.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Constructor for the xyz coordinate data.
 *  \return See brief. */
static const struct const_Multiarray_d* constructor_xyz ();

/** \brief Constructor for the normal vector data.
 *  \return See brief. */
/// \todo make static
const struct const_Multiarray_d* constructor_normals ();

/// \brief Set the members of the input real \ref mutable_Flux_T container to zero.
static void set_to_zero_Flux
	(struct mutable_Flux*const m_flux,    ///< \ref mutable_Flux_T.
	 const struct Flux_Input*const flux_i ///< \ref Flux_Input_T.
	);

static void constructor_Flux_Input_data_members (struct Flux_Input*const flux_i, const struct Simulation*const sim)
{
	flux_i->s   = constructor_s(0,sim);
	flux_i->g   = constructor_g(0,sim);
	flux_i->xyz = constructor_xyz();
}

static void destructor_Flux_Input_data_members (struct Flux_Input*const flux_i)
{
	destructor_const_Multiarray_d(flux_i->s);
	destructor_const_Multiarray_d(flux_i->g);
	destructor_const_Multiarray_d(flux_i->xyz);
}

static void constructor_Flux_Input_c_data_members (struct Flux_Input_c*const flux_c_i, struct Flux_Input*const flux_i)
{
	flux_c_i->s   = constructor_copy_const_Multiarray_c_Multiarray_d(flux_i->s);
	flux_c_i->g   = constructor_copy_const_Multiarray_c_Multiarray_d(flux_i->g);
	flux_c_i->xyz = constructor_copy_const_Multiarray_d(flux_i->xyz);
}

static void destructor_Flux_Input_c_data_members (struct Flux_Input_c*const flux_c_i)
{
	destructor_const_Multiarray_c(flux_c_i->s);
	destructor_const_Multiarray_c(flux_c_i->g);
	destructor_const_Multiarray_d(flux_c_i->xyz);
}

static struct Flux* constructor_Flux_cmplx_step (struct Flux_Input_c*const flux_i_c, struct Flux_Input*const flux_i)
{
	struct Flux* flux = constructor_Flux(flux_i); // returned
	struct mutable_Flux*const m_flux = (struct mutable_Flux*) flux;
	set_to_zero_Flux(m_flux,flux_i);

	const bool* compute_member = flux_i->compute_member;

	struct Multiarray_c* s = (struct Multiarray_c*) flux_i_c->s;
	const ptrdiff_t n_n  = s->extents[0],
	                n_eq = flux_i->n_eq,
	                n_vr = flux_i->n_var;

	for (int vr = 0; vr < n_vr; ++vr) {
		double complex*const data_s = get_col_Multiarray_c(vr,s);
		add_to_c(data_s,CX_STEP*I,n_n);
		struct Flux_c* flux_c = constructor_Flux_c(flux_i_c); // destructed

		// f[DIM,NEQ]
		assert(compute_member[0] == true);
		if (vr == 0) {
			for (int eq = 0; eq < n_eq; ++eq) {
			for (int d = 0; d < DIM; ++d) {
				const ptrdiff_t ind_c = d+DIM*(eq),
				                ind_r = d+DIM*(eq);

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,flux_c->f);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_flux->f);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = creal(data_c[n]);
			}}
		}

		// df_ds[DIM,NEQ,NVAR]
		if (compute_member[1]) {
			for (int eq = 0; eq < n_eq; ++eq) {
			for (int d = 0; d < DIM; ++d) {
				const ptrdiff_t ind_c = d+DIM*(eq),
				                ind_r = d+DIM*(eq+n_eq*(vr));

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,flux_c->f);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_flux->df_ds);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}}
		}

		// d2f_ds2[DIM,NEQ,NVAR,NVAR]
		if (compute_member[3]) {
			for (int vr2 = 0; vr2 < n_vr; ++vr2) {
			for (int eq = 0; eq < n_eq; ++eq) {
			for (int d = 0; d < DIM; ++d) {
				const ptrdiff_t ind_c = d+DIM*(eq+n_eq*(vr2)),
				                ind_r = d+DIM*(eq+n_eq*(vr2+n_vr*(vr)));

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,flux_c->df_ds);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_flux->d2f_ds2);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}}}
		}

		destructor_Flux_c(flux_c);
		add_to_c(data_s,-CX_STEP*I,n_n);
	}

	struct Multiarray_c*const g = (struct Multiarray_c*) flux_i_c->g;

	for (int dx = 0; dx < DIM; ++dx) {
	for (int vr = 0; vr < n_vr; ++vr) {
		double complex*const data_g = get_col_Multiarray_c(vr+n_vr*(dx),g);
		add_to_c(data_g,CX_STEP*I,n_n);
		struct Flux_c* flux_c = constructor_Flux_c(flux_i_c); // destructed

		// df_dg[DIM,NEQ,NVAR,DIM]
		if (compute_member[2]) {
			for (int eq = 0; eq < n_eq; ++eq) {
			for (int d = 0; d < DIM; ++d) {
				const ptrdiff_t ind_c = d+DIM*(eq),
				                ind_r = d+DIM*(eq+n_eq*(vr+n_vr*(dx)));

				const double complex*const data_c = get_col_const_Multiarray_c(ind_c,flux_c->f);
				double*const data_r               = get_col_Multiarray_d(ind_r,m_flux->df_dg);
				for (int n = 0; n < n_n; ++n)
					data_r[n] = cimag(data_c[n])/CX_STEP;
			}}
		}

		assert(!compute_member[4]); // Add support.
		assert(!compute_member[5]); // Add support.

		destructor_Flux_c(flux_c);
		add_to_c(data_g,-CX_STEP*I,n_n);
	}}

	return flux;
}

// Level 1 ********************************************************************************************************** //

static const struct const_Multiarray_d* constructor_s (const int index, const struct Simulation*const sim)
{
	const double* data = NULL;
	switch (DIM) {
	case 1:
		switch (index) {
		case 0: {
			static const double s[] = {  1.01,  1.02,  1.03,
			                             2.01, -2.02,  2.03,
			                             9.01,  9.02,  9.03, };
			data = s;
			break;
		} case 1: {
			static const double s[] = {  1.11,  1.12,  1.13,
			                            -2.11,  2.12,  2.53,
			                             9.11,  9.12,  9.13, };
			data = s;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	case 2:
		switch (index) {
		case 0: {
			static const double s[] = {  1.01,  1.02,  1.03,
			                             2.01, -2.22,  2.03,
			                             2.04, -2.25,  2.06,
			                             9.01,  9.02,  9.03, };
			data = s;
			break;
		} case 1: {
			static const double s[] = {  1.11,  1.12,  1.13,
			                            -2.11,  2.12,  2.53,
			                            -2.14,  2.15,  2.56,
			                             9.11,  9.12,  9.13, };
			data = s;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	case 3:
		switch (index) {
		case 0: {
			static const double s[] = {  1.01,  1.02,  1.03,
			                             2.41, -2.02,  1.92,
			                             2.44, -2.05,  1.95,
			                             2.47, -2.08,  1.98,
			                             9.01,  9.02,  9.02, };
			data = s;
			break;
		} case 1: {
			static const double s[] = {  1.11,  1.12,  1.13,
			                            -2.11,  2.12,  2.13,
			                            -2.14,  2.15,  2.16,
			                            -2.17,  2.18,  2.19,
			                             9.11,  9.12,  9.13, };
			data = s;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",DIM);
		break;
	}

	struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;

	const ptrdiff_t* extents = (ptrdiff_t[]) { 3, test_case->n_var, };
	return constructor_move_const_Multiarray_d_d('C',2,extents,false,data);
}

static const struct const_Multiarray_d* constructor_g (const int index, const struct Simulation*const sim)
{
	const double* data = NULL;
	switch (DIM) {
	case 1:
		switch (index) {
		case 0: {
			static const double g[] = { -0.0855,  -0.9289,   0.2373,
			                            -0.2625,   0.7303,   0.4588,
			                            -0.0292,  -0.5785,  -0.5468, };
			data = g;
			break;
		} case 1: {
			static const double g[] = { -0.5211,   0.6791,  -0.0377,
			                            -0.2316,   0.3955,  -0.8852,
			                            -0.6241,  -0.9880,   0.7962, };
			data = g;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	case 2:
		switch (index) {
		case 0: {
			static const double g[] = { -0.0855,  -0.9289,   0.2373,
			                            -0.2625,   0.7303,   0.4588,
			                            -0.8010,  -0.4886,  -0.9631,
			                            -0.0292,  -0.5785,  -0.5468,

			                            -0.4018,   0.1839,  -0.9027,
			                             0.0760,   0.2400,   0.9448,
			                            -0.2399,  -0.4173,  -0.4909,
			                             0.1233,  -0.0497,   0.4893, };
			data = g;
			break;
		} case 1: {
			static const double g[] = { -0.5211,   0.6791,  -0.0377,
			                            -0.2316,   0.3955,  -0.8852,
			                             0.4889,  -0.3674,  -0.9133,
			                            -0.6241,  -0.9880,   0.7962,

			                            -0.3377,   0.7803,   0.0965,
			                            -0.9001,   0.3897,   0.1320,
			                            -0.3692,   0.2417,  -0.9421,
			                            -0.1112,  -0.4039,   0.9561, };
			data = g;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	case 3:
		switch (index) {
		case 0: {
			static const double g[] = { -0.0911,  -0.6444,  -0.2089,
			                             0.5762,   0.6476,   0.7093,
			                             0.6834,  -0.6790,  -0.2362,
			                            -0.5466,  -0.6358,   0.1194,
			                             0.4257,  -0.9452,  -0.6073,

			                             0.7184,  -0.6110,   0.1537,
			                             0.9686,   0.7788,  -0.2810,
			                             0.5313,   0.4235,  -0.4401,
			                            -0.3251,   0.0908,  -0.5271,
			                             0.1056,  -0.2665,  -0.4574,

			                            -0.4899,   0.4711,   0.5216,
			                            -0.1679,  -0.0596,   0.0967,
			                            -0.9787,   0.6820,   0.8181,
			                             0.7127,  -0.0424,  -0.8175,
			                            -0.5005,   0.0714,   0.7224, };
			data = g;
			break;
		} case 1: {
			static const double g[] = { -0.4501,  -0.6620,   0.6135,
			                            -0.4587,   0.4162,  -0.5822,
			                             0.6619,  -0.8419,   0.5407,
			                            -0.7703,   0.8329,  -0.8699,
			                            -0.3502,  -0.2564,  -0.2648,

			                            -0.8754,   0.2407,  -0.0680,
			                            -0.5181,   0.6761,   0.2548,
			                            -0.9436,   0.2891,   0.2240,
			                             0.6377,  -0.6718,  -0.6678,
			                            -0.9577,  -0.6951,   0.8444,

			                            -0.1499,   0.8003,   0.1332,
			                            -0.6596,   0.4538,  -0.1734,
			                             0.5186,   0.4324,  -0.3909,
			                             0.9730,   0.8253,   0.8314,
			                            -0.6490,   0.0835,  -0.8034, };
			data = g;
			break;
		} default:
			EXIT_ERROR("Unsupported: %d\n",index);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",DIM);
		break;
	}

	struct Test_Case* test_case = (struct Test_Case*) sim->test_case_rc->tc;

	const ptrdiff_t* extents = (ptrdiff_t[]) { 3, test_case->n_var, DIM };
	return constructor_move_const_Multiarray_d_d('C',3,extents,false,data);
}

static const struct const_Multiarray_d* constructor_xyz ()
{
	static const double xyz[] = {  1.8147,  1.8568,  0.9360,
	                               0.2785,  0.8559,  1.5453,
	                               0.1143, -0.0037,  0.0751, };
	const double* data = xyz;

	const ptrdiff_t* extents = (ptrdiff_t[]) { 3, DIM };
	return constructor_move_const_Multiarray_d_d('C',2,extents,false,data);
}

const struct const_Multiarray_d* constructor_normals ()
{
	const double* data = NULL;
	switch (DIM) {
	case 1: {
		static const double n[] = {  0.6651,
		                             0.8190,
		                             0.2449, };
		data = n;
		break;
	} case 2: {
		static const double n[] = {  0.6651,  0.7395,
		                             0.8190,  0.5670,
		                             0.2449, -0.4809, };
		data = n;
		break;
	} case 3: {
		static const double n[] = {  0.6651,  0.7395,  0.1037,
		                             0.8190,  0.5670, -0.0875,
		                             0.8190,  0.5670, -0.0875, };
		data = n;
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",DIM);
		break;
	}

	const ptrdiff_t* extents = (ptrdiff_t[]) { 3, DIM };
	return constructor_move_const_Multiarray_d_d('R',2,extents,false,data);
}

static void set_to_zero_Flux (struct mutable_Flux*const m_flux, const struct Flux_Input*const flux_i)
{
	const bool*const c_m = flux_i->compute_member;

	if (c_m[0])
		set_to_value_Multiarray_d(m_flux->f,0.0);
	if (c_m[1])
		set_to_value_Multiarray_d(m_flux->df_ds,0.0);
	if (c_m[2])
		set_to_value_Multiarray_d(m_flux->df_dg,0.0);
	if (c_m[3])
		set_to_value_Multiarray_d(m_flux->d2f_ds2,0.0);
	if (c_m[4])
		EXIT_ADD_SUPPORT;
	if (c_m[5])
		EXIT_ADD_SUPPORT;
}
