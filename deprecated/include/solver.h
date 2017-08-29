// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_h__INCLUDED
#define DPG__solver_h__INCLUDED

#include <stdbool.h>
#include "petscmat.h"

#include "S_DB.h" // For constructor_Simulation (ToBeDeleted)
#include "S_VOLUME.h"

struct S_solver_info {
	/*
	 *	Purpose:
	 *		Container used to pass around solver-related information.
	 *
	 *	Notation:
	 *		display    : display information while solving.
	 *		output     : output the solution (for visualization in paraview) while solving.
	 *		adapt      : update the finite element space (through hp adaptation) while solving.
	 *		positivity : enforce positivity of the solution based on physical constraints.
	 *		symmetric  : use a symmetric implicit solver procedure.
	 *		steady     : converge to a steady-state solution (as opposed to an unsteady solution).
	 *		linear     : indication that the PDE is linear, meaning that steady solutions can be computed in a single
	 *		             Newton step.
	 *		create_RHS : create petsc Vec's for x and b (in Ax = b). Enabled by default. Disabled for linearization
	 *		             testing.
	 *		compute_V  : Specify whether (V)OLUME terms should be included in the computation. Enabled by default.
	 *		             Disabled for partial linearization testing.
	 *		compute_F  : Specify whether FACE terms should be included in the computation. Enabled by default.
	 *		              Disabled for partial linearization testing.
	 *
	 *		imex_type : (im)plicit-(ex)plicit (type) (i.e. whether using an implicit or explicit solver).
	 *		method    : method used for the computation of the solution.
	 *
	 *		dof       : (d)egrees (o)f (f)reedom in the global system.
	 *		nnz       : (n)umber of (n)on-(z)ero entries in each row of the global system matrix.
	 *		A, x, b   : Petsc containers for global system solve (Ax = b).
	 *
	 *		compute_all      : Specify whether contributions from all VOLUMEs should be computed. Disabled by default.
	 *		                   Enabled for testing equivalence of real and complex functions.
	 *		VOLUME_perturbed : (VOLUME) which has been (perturbed) with the complex step.
	 */

	bool display, output, adapt, positivity, symmetric, steady, linear, create_RHS, compute_V, compute_F;
	char imex_type;
	unsigned int method;

	// Petsc related variables (Only used for linearization)
	unsigned int dof;
	PetscInt *nnz;

	Mat A;
	Vec x,
	    b;

	// Complex-step verification variables (Only used for testing)
	bool compute_all;
	struct S_VOLUME *VOLUME_perturbed;
};

struct S_LHS_info {
	size_t IndA[2], Nn[2];
	double* LHS;

	const struct S_VOLUME* VOLUME[2];

	InsertMode addv;
};

extern struct S_solver_info constructor_solver_info (const bool display, const bool output, const bool adapt,
                                                     const char imex_type, const unsigned int method);
extern void set_global_indices       (struct S_solver_info*const solver_info);
extern void initialize_petsc_structs (struct S_solver_info*const solver_info);
extern void assemble_petsc_structs   (struct S_solver_info*const solver_info);
extern void destroy_petsc_structs    (struct S_solver_info*const solver_info);

extern struct S_LHS_info constructor_LHS_info (double*const LHS, const struct S_VOLUME*const V0,
                                               const struct S_VOLUME*const V1, const InsertMode addv,
                                               const bool correct_symm);

extern void compute_GradW (const struct S_solver_info*const solver_info, const char stage);
extern void compute_RLHS  (const struct S_solver_info*const solver_info);
extern void fill_PetscMat (const struct S_solver_info*const solver_info, const struct S_LHS_info*const LHS_info);

#endif // DPG__solver_h__INCLUDED
