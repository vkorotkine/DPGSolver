// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_h__INCLUDED
#define DPG__solver_h__INCLUDED

#include <stdbool.h>
#include "petscmat.h"

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
	 *
	 *		imex_type : (im)plicit-(ex)plicit (type) (i.e. whether using an implicit or explicit solver).
	 *		method    : method used for the computation of the solution.
	 */

	bool display, output, adapt, positivity, symmetric, steady, linear;
	char imex_type;
	unsigned int method;

	// Petsc related parameters
	unsigned int dof;
	PetscInt *nnz;

	Mat A;
	Vec x,
	    b;
};

struct S_LHS_info {
	size_t IndA[2], Nn[2];
	const double* LHS;

	InsertMode addv;
};

extern struct S_solver_info constructor_solver_info (const bool display, const bool output, const bool adapt,
                                                     const char imex_type, const unsigned int method);
extern void initialize_petsc_structs (struct S_solver_info*const solver_info);

extern struct S_LHS_info constructor_LHS_info (const double*const LHS, const struct S_VOLUME*const V0,
                                               const struct S_VOLUME*const V1, const InsertMode addv);

extern void compute_RLHS  (const struct S_solver_info*const solver_info);
extern void fill_PetscMat (const struct S_solver_info*const solver_info, const struct S_LHS_info*const LHS_info);

#endif // DPG__solver_h__INCLUDED
