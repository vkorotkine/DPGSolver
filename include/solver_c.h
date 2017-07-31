// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__solver_c_h__INCLUDED
#define DPG__solver_c_h__INCLUDED

#include "containers.h"
#include "containers_c.h"
#include "solver.h"

struct Volume_solver_c {
	// Moved variables
	const struct const_Multiarray_d*const XYZ,
	                               *const detJV_vI,
	                               *const C_vI;

	// Function specific allocated variables
	struct Multiarray_c*const What,  // Solution coefficients
	                   *const QhatV, // Local weak gradient coefficients
	                   *const Qhat,  // Weak gradient coefficients
	                   *const RHS;   // Residual

	// Additional variables
	struct Volume_solver_c* next;
};

struct Context_solver_c {
	const unsigned int d,
	                   n_var;

	struct Volume_solver_c*const volume_head; // Pointer to first Volume_solver_c
};

struct Context_solver_c constructor_Context_solver_c
	(const struct Simulation*const simulation, const struct S_VOLUME*const VOLUME_head);

extern void compute_GradW_c (const struct S_solver_info*const solver_info, const char stage);

#endif // DPG__solver_c_h__INCLUDED
