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

#ifndef DPG__solution_navier_stokes_h__INCLUDED
#define DPG__solution_navier_stokes_h__INCLUDED
/** \file
 *  \brief Provides real functions relating to the Navier-Stokes solutions.
 */

#include "def_templates_type_d.h"
#include "def_templates_solution_navier_stokes.h"
#include "def_templates_multiarray.h"
#include "def_templates_test_case.h"
#include "solution_navier_stokes_T.h"
#include "undef_templates_type.h"
#include "undef_templates_solution_navier_stokes.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_test_case.h"

/** \brief Return the value of the coefficient of thermal conductivity for the Navier-Stokes equations, assuming
 *         constant specific heat at constant pressure (eq. (1.59), \cite Toro2009).
 *  \return See brief. */
double compute_kappa_const_cp
	(const double mu, ///< The viscosity.
	 const double Cp, ///< The heat capacity at constant pressure.
	 const double Pr  ///< The Prandtl number.
	);

/** \brief Return the value of the heat capacity at constant pressure for an ideal gas (eq. (1.46), \cite Toro2009).
 *  \return See brief. */
double compute_cp_ideal_gas
	(const double r_s ///< The specific gas constant.
	);

/// \brief Compute the viscosity.
void compute_viscosity
	(struct Multiarray_d* mu,               ///< The container to hold the viscosity data.
	 const struct const_Multiarray_d* vars, ///< The container of Euler variables.
	 const char var_type                    ///< The type of the variables.
	);

/** \brief Get the normal flux for the Energy equation in cases where the value is constant along an entire boundary.
 *  \return See brief. */
double get_normal_flux_Energy
	( );

#endif // DPG__solution_navier_stokes_h__INCLUDED
