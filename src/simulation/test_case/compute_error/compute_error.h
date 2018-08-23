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

#ifndef DPG__compute_error_h__INCLUDED
#define DPG__compute_error_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for error computation and output.
 */

#include "compute_error_advection.h"
#include "compute_error_diffusion.h"
#include "compute_error_euler.h"
#include "compute_error_navier_stokes.h"

struct Error_CE;
struct Simulation;
struct Solver_Volume;
struct Vector_d;
struct const_Multiarray_d;

#include <stddef.h>
#include "definitions_alloc.h"

/** \brief Function pointer to error computing functions (which construct \ref Error_CE containers).
 *  \param sim \ref Simulation.
 */
typedef struct Error_CE* (*constructor_Error_CE_fptr)
	(const struct Simulation* sim
	);

/// Container holding information relating to the computational element errors.
struct Error_CE {
	const ptrdiff_t dof;        ///< The number of degrees of freedom.
	const double domain_volume; ///< The volume of the domain.

	const struct const_Vector_d* sol_err;        ///< The solution error of the desired error type.
	const struct const_Vector_i* expected_order; ///< The expected convergence order of each variable.

	const char* header_spec; ///< The header for the file specific to the variables for which the error was computed.
};

/// Container holding information used to construct the \ref Error_CE container.
struct Error_CE_Helper {
	const int n_out;      ///< The number of output error variables.
	const int error_type; ///< The type of error. Options: see \ref definitions_error.h
	int domain_order;     ///< Polynomial order of computational elements in the domain.
	double domain_volume; ///< \ref Error_CE::domain_volume.

	const char* header_spec; ///< The header for the file specific to the variables for which the error was computed.
	struct Vector_d* sol_err;        ///< \ref Error_CE::sol_err.

	struct Solution_Container* sol_cont; ///< \ref Solution_Container_T.

	struct Solver_Volume* s_vol[2];   ///< Pointers to the computed and exact \ref Solver_Volume_T\*s.
	const struct Solver_Face* s_face; ///< Pointer to the \ref Solver_Face_T (if applicable).
};

/// Container holding information relating to the solution used to compute the error.
struct Error_CE_Data {
	struct Multiarray_d* sol[2]; ///< Pointers to the computed and exact solution data.
	struct Multiarray_d* rhs[2]; ///< Pointers to the computed and exact rhs data.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_Error_CE_fptr returning NULL;
 *  \return See brief. */
struct Error_CE* constructor_Error_CE_NULL
	(const struct Simulation* sim ///< Defined for \ref constructor_Error_CE_fptr.
	 );

/** \brief Output the error of the solution.
 *
 *  The error, for a generic variable \f$v\f$, is computed as:
 *  \f[
 *  ||v||_{L^2(\Omega)} =
 *  \left(
 *  \frac{\int_{\Omega} (v-v_{\text{exact}})^2 d \Omega}{\int_{\Omega} d \Omega}
 *  \right)^{\frac{1}{2}}
 *  \f]
 */
void output_error
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output the error of the functionals.
void output_error_functionals
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Error_CE container.
 *  \return See brief. */
struct Error_CE* constructor_Error_CE
	(struct Error_CE_Helper* e_ce_h, ///< \ref Error_CE_Helper.
	 const struct Simulation* sim    ///< \ref Simulation.
	);

/** \brief Constructor for a \ref Error_CE_Helper container.
 *  \return See brief. */
struct Error_CE_Helper* constructor_Error_CE_Helper
	(const struct Simulation* sim, ///< \ref Simulation.
	 const int n_out               ///< The number of output error variables.
	);

/// \brief Destructor for a \ref Error_CE_Helper container.
void destructor_Error_CE_Helper
	(struct Error_CE_Helper* e_ce_h ///< Standard.
	);

/** \brief Constructor for a \ref Error_CE_Data container.
 *  \return See brief. */
struct Error_CE_Data* constructor_Error_CE_Data
	(struct Error_CE_Helper* e_ce_h, ///< \ref Error_CE_Helper.
	 const struct Simulation* sim    ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Error_CE_Data container.
void destructor_Error_CE_Data
	(struct Error_CE_Data* e_ce_d ///< Standard.
	);

/// \brief Increment \ref Error_CE_Helper::sol_err.
void increment_sol_L2
	(struct Error_CE_Helper* e_ce_h, ///< \ref Error_CE_Helper.
	 struct Error_CE_Data* e_ce_d    ///< \ref Error_CE_Data.
	);

/// \brief Increment \ref Error_CE_Helper::sol_err for the input face variables.
void increment_sol_integrated_face
	(struct Error_CE_Helper*const e_ce_h, ///< \ref Error_CE_Helper.
	 struct Error_CE_Data*const e_ce_d    ///< \ref Error_CE_Data.
	);

/// \brief Increment \ref Error_CE_Helper::sol_err on a face.
void increment_sol_face_L2
	(struct Error_CE_Helper*const e_ce_h, ///< \ref Error_CE_Helper.
	 struct Error_CE_Data*const e_ce_d    ///< \ref Error_CE_Data.
	);

/// \brief Update \ref Error_CE_Helper::domain_order.
void update_domain_order
	(struct Error_CE_Helper* e_ce_h ///< \ref Error_CE_Helper.
	);

/** \brief Return a statically allocated `char*` holding the name of the error file.
 *  \return See brief. */
const char* compute_error_file_name
	(const int error_type,             ///< The number representing the type of error being output.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Correct existing (add if missing) 'ml' and 'p' parameters of the file name if input values are valid.
void correct_file_name_ml_p
	(const int ml,        ///< The mesh level.
	 const int p,         ///< The polynomial order.
	 char*const file_name ///< The input file name.
	);

/// \brief Add the computed and exact specified variable to the solution data Multiarrays.
void add_rhs_Error_CE_Data
	(struct Error_CE_Data*const e_ce_d, ///< \ref Error_CE_Data.
	 const struct Simulation*const sim  ///< \ref Simulation.
	);

#endif // DPG__compute_error_h__INCLUDED
