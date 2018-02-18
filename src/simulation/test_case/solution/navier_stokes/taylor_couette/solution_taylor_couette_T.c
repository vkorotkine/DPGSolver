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
#include <math.h>
#include <string.h>

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_intrusive.h"
#include "definitions_math.h"
#include "definitions_physics.h"
#include "definitions_tol.h"


#include "def_templates_solution.h"
#include "def_templates_solution_navier_stokes.h"

#include "def_templates_multiarray.h"

#include "def_templates_math_functions.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Return a \ref Multiarray_T\* container holding the solution values at the input coordinates.
 *  \return See brief.
 *
 *  \note Only the velocity and temperature components are exact.
 */
static struct Multiarray_T* constructor_sol_taylor_couette
	(const struct const_Multiarray_R*const xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation*const sim          ///< \ref Simulation.
	);

/** \brief Return a \ref Multiarray_T\* container holding the solution gradient values at the input coordinates.
 *  \return See brief.
 *
 *  \note Only the velocity and temperature components are exact.
 */
static struct Multiarray_T* constructor_grad_taylor_couette
	(const struct const_Multiarray_R*const xyz, ///< xyz coordinates at which to evaluate the solution.
	 const struct Simulation*const sim          ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void set_sol_taylor_couette_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(sol_cont.node_kind != 'r');
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* sol = constructor_sol_taylor_couette(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_sol_T(&sol_cont,sol,sim);
	destructor_Multiarray_T(sol);
}

void set_grad_taylor_couette_T (const struct Simulation* sim, struct Solution_Container_T sol_cont)
{
	assert(sol_cont.node_kind != 's');
	const struct const_Multiarray_R* xyz = constructor_xyz_sol_T(sim,&sol_cont); // destructed
	struct Multiarray_T* grad = constructor_grad_taylor_couette(xyz,sim); // destructed
	destructor_const_Multiarray_R(xyz);

	update_Solution_Container_grad_T(&sol_cont,grad,sim);
	destructor_Multiarray_T(grad);
}

const struct const_Multiarray_T* constructor_const_sol_taylor_couette_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	return (struct const_Multiarray_T*) constructor_sol_taylor_couette(xyz,sim);
}

const struct const_Multiarray_T* constructor_const_grad_taylor_couette_T
	(const struct const_Multiarray_R* xyz, const struct Simulation* sim)
{
	return (struct const_Multiarray_T*) constructor_grad_taylor_couette(xyz,sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for solution data relating to 't'aylor-'c'ouette.
struct Sol_Data__tc {
	// Read parameters
	Real r_i,    ///< The 'r'adius of the 'i'nner cylinder.
	     r_o,    ///< The 'r'adius of the 'o'uter cylinder.

	     Pr,     ///< The Prandtl number.
	     mu,     ///< The constant viscosity value.
	     omega,  ///< The angular velocity of the rotating cylinder.
	     dTdr_o, ///< The gradient of the temperature in the radial direction at the 'o'uter cylinder boundary.
	     t_i,    ///< The temperature of the 'i'nner cylinder.
	     p_i,    ///< The pressure of the 'i'nner cylinder.
	     r_s;    ///< The specific gas constant.

	// Additional parameters
	Real kappa; ///< The coefficient of thermal conductivity.
};

/** \brief Return the statically allocated \ref Sol_Data__tc container.
 *  \return See brief. */
static struct Sol_Data__tc get_sol_data
	(const struct Simulation*const sim ///< \ref Simulation.
	);

static struct Multiarray_T* constructor_sol_taylor_couette
	(const struct const_Multiarray_R*const xyz, const struct Simulation*const sim)
{
	assert(DIM == 2);
	const struct Sol_Data__tc sol_data = get_sol_data(sim);

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	const Real*const x = get_col_const_Multiarray_R(0,xyz),
	          *const y = get_col_const_Multiarray_R(1,xyz);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	struct Multiarray_T* sol = constructor_empty_Multiarray_T('C',2,(ptrdiff_t[]){n_n,n_var}); // returned

	Type*const rho_ptr = get_col_Multiarray_T(0,sol),
	    *const u_ptr   = get_col_Multiarray_T(1,sol),
	    *const v_ptr   = get_col_Multiarray_T(2,sol),
	    *const p_ptr   = get_col_Multiarray_T(n_var-1,sol);

	const Real r_i    = sol_data.r_i,
	           r_o    = sol_data.r_o,
	           omega  = sol_data.omega,
	           dTdr_o = sol_data.dTdr_o,
	           t_i    = sol_data.t_i,
	           p_i    = sol_data.p_i,
	           mu     = sol_data.mu,
	           kappa  = sol_data.kappa,
	           r_s    = sol_data.r_s,

	           c      = omega/(1.0/(r_i*r_i)-1.0/(r_o*r_o));

	for (int i = 0; i < n_n; ++i) {
		const Real r  = sqrt(x[i]*x[i]+y[i]*y[i]),
		           th = atan2(y[i],x[i]),
		           Vt = c*(1.0/r-r/(r_o*r_o)),
		           t  = t_i - 2.0*c*c/(r_o*r_o)*mu/kappa*log(r/r_i) - c*c*mu/kappa*(1.0/(r*r)-1.0/(r_i*r_i))
		              + dTdr_o*r_o*log(r/r_i),
		           p  = p_i;

		u_ptr[i]   = -sin(th)*Vt;
		v_ptr[i]   =  cos(th)*Vt;

		/** \note The pressure is nearly uniform (p.8, \cite Illingworth1950), hence the pressure is set to that of
		 *        the internal cylinder and the density is computed using the equation of state. */
		p_ptr[i]   = p;
		rho_ptr[i] = p/(r_s*t);
	}
	convert_variables_T(sol,'p','c');

	return sol;
}

static struct Multiarray_T* constructor_grad_taylor_couette
	(const struct const_Multiarray_R*const xyz, const struct Simulation*const sim)
{
	assert(DIM == 2);
	static const int dmax = 2;

	const struct Sol_Data__tc sol_data = get_sol_data(sim);

	// Compute the solution
	const ptrdiff_t n_n = xyz->extents[0];
	assert(DIM == xyz->extents[1]);

	const Real*const x = get_col_const_Multiarray_R(0,xyz),
	          *const y = get_col_const_Multiarray_R(1,xyz);

	struct Test_Case_T* test_case = (struct Test_Case_T*)sim->test_case_rc->tc;
	const int n_var = test_case->n_var;

	const struct const_Multiarray_T*const sol = (struct const_Multiarray_T*)
		constructor_sol_taylor_couette(xyz,sim); // destructed
	const Type*const rho_ptr = get_col_const_Multiarray_T(0,sol),
	          *const u_ptr   = get_col_const_Multiarray_T(1,sol),
	          *const v_ptr   = get_col_const_Multiarray_T(2,sol),
	          *const p_ptr   = get_col_const_Multiarray_T(n_var-1,sol);

	struct Multiarray_T*const grad = constructor_empty_Multiarray_T('C',3,(ptrdiff_t[]){n_n,n_var,DIM}); // returned
	Type*const g_rho_ptr[] = ARRAY_DIM( get_col_Multiarray_T(0*n_var+0,grad), get_col_Multiarray_T(1*n_var+0,grad), NULL ),
	    *const g_u_ptr[]   = ARRAY_DIM( get_col_Multiarray_T(0*n_var+1,grad), get_col_Multiarray_T(1*n_var+1,grad), NULL ),
	    *const g_v_ptr[]   = ARRAY_DIM( get_col_Multiarray_T(0*n_var+2,grad), get_col_Multiarray_T(1*n_var+2,grad), NULL ),
	    *const g_p_ptr[]   = ARRAY_DIM( get_col_Multiarray_T(0*n_var+3,grad), get_col_Multiarray_T(1*n_var+3,grad), NULL );

	const Real r_i    = sol_data.r_i,
	           r_o    = sol_data.r_o,
	           omega  = sol_data.omega,
	           dTdr_o = sol_data.dTdr_o,
	           mu     = sol_data.mu,
	           kappa  = sol_data.kappa,
	           r_s    = sol_data.r_s,

	           c      = omega/(1.0/(r_i*r_i)-1.0/(r_o*r_o));

	for (int i = 0; i < n_n; ++i) {
		const Real r  = sqrt(x[i]*x[i]+y[i]*y[i]),
		           th = atan2(y[i],x[i]);
		const Type Vt = u_ptr[i]/(-sin(th)),
		           t  = p_ptr[i]/(r_s*rho_ptr[i]),
			     p  = p_ptr[i];
		assert(equal_T(Vt,v_ptr[i]/(cos(th)),EPS));

		const Real dr_dX[]  = { x[i]/r, y[i]/r, },
		           dth_dX[] = { -y[i]/(r*r), x[i]/(r*r), };

		const Real dVt_dr = c*(-1.0/(r*r)-1.0/(r_o*r_o)),
		           dt_dr  = -c*c*mu/kappa*(2.0/r*(1.0/(r_o*r_o)-1.0/(r*r)))+dTdr_o*r_o/r;

		/** \note The values for the gradients of the density and pressure are not exact, but must be consistent
		 *        with the equation of state, such that the correct temperature gradient is obtained. */
		const Real dp_dX[] = { 0.0, 0.0, };

		for (int d = 0; d < dmax; ++d) {
			g_rho_ptr[d][i] = 1.0/(r_s*t)*dp_dX[d] - p/(r_s*t*t)*dt_dr*dr_dX[d];
			g_u_ptr[d][i]   = -(  cos(th)*dth_dX[d]*Vt + sin(th)*dVt_dr*dr_dX[d] );
			g_v_ptr[d][i]   =  ( -sin(th)*dth_dX[d]*Vt + cos(th)*dVt_dr*dr_dX[d] );
			g_p_ptr[d][i]   = dp_dX[d];
		}
	}
	convert_variables_gradients_T(grad,sol,'p','c');
	destructor_const_Multiarray_T(sol);

	return grad;
}

// Level 1 ********************************************************************************************************** //

/// \brief Read the required solution data into \ref Sol_Data__tc.
static void read_data_taylor_couette
	(struct Sol_Data__tc*const sol_data ///< \ref Sol_Data__tc.
	);

/// \brief Set the remaining required solution data of \ref Sol_Data__tc based on the read values.
static void set_data_taylor_couette
	(struct Sol_Data__tc*const sol_data ///< \ref Sol_Data__tc.
	);

static struct Sol_Data__tc get_sol_data (const struct Simulation*const sim)
{
	static bool need_input = true;

	static struct Sol_Data__tc sol_data;
	if (need_input) {
		need_input = false;
		read_data_taylor_couette(&sol_data);
		set_data_taylor_couette(&sol_data);

		if (!strstr(sim->geom_spec,"diabatic_o"))
			EXIT_ERROR("Unsupported: %s\n",sim->geom_spec);
	}

	return sol_data;
}

// Level 2 ********************************************************************************************************** //

static void read_data_taylor_couette (struct Sol_Data__tc*const sol_data)
{
	const int count_to_find = 11;
	int count_found = 0;

	FILE* input_file = NULL;
	char line[STRLEN_MAX];

	input_file = fopen_input('g',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_string_count_c_style_d("r_i",&count_found,line,&sol_data->r_i);
		read_skip_string_count_c_style_d("r_o",&count_found,line,&sol_data->r_o);
	}
	fclose(input_file);

	int viscosity_type = VISCOSITY_INVALID;
	int diabatic_flux_type = VISCOUS_BC_INVALID;

	input_file = fopen_input('s',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_convert_i(line,"viscosity_type",&viscosity_type,&count_found);
		read_skip_convert_i(line,"diabatic_flux_type",&diabatic_flux_type,&count_found);

		read_skip_string_count_d("Pr",&count_found,line,&sol_data->Pr);
		read_skip_string_count_d("mu",&count_found,line,&sol_data->mu);

		read_skip_string_count_d("omega",&count_found,line,&sol_data->omega);
		read_skip_string_count_d("dTdn" ,&count_found,line,&sol_data->dTdr_o);
		read_skip_string_count_d("t_b",  &count_found,line,&sol_data->t_i);
		read_skip_string_count_d("p_b",  &count_found,line,&sol_data->p_i);
		read_skip_string_count_d("r_s",  &count_found,line,&sol_data->r_s);
	}
	fclose(input_file);
	assert(count_found == count_to_find);

	assert(viscosity_type == VISCOSITY_CONSTANT);
	assert(((diabatic_flux_type == DIABATIC_FLUX_CONSTANT_ZERO) && sol_data->dTdr_o == 0.0) ||
	       ((diabatic_flux_type == DIABATIC_FLUX_CONSTANT)      && sol_data->dTdr_o != 0.0));
}

static void set_data_taylor_couette (struct Sol_Data__tc*const sol_data)
{
	const Real Cp = compute_cp_ideal_gas(sol_data->r_s);
	sol_data->kappa = compute_kappa_const_cp(sol_data->mu,Cp,sol_data->Pr);
}
