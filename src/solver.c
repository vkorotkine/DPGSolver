// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "solver.h"
#include <stdbool.h>
#include "petscmat.h"
#include "S_VOLUME.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_FACE.h"

#include "compute_VOLUME_RLHS_HDG.h"
#include "compute_FACE_RLHS_HDG.h"

#include "compute_GradW_DG.h"
#include "compute_VOLUME_RLHS_DG.h"
#include "compute_FACE_RLHS_DG.h"
#include "solver_symmetric_functions.h"

/*
 *	Purpose:
 *		Provide solver related functions.
 *
 *	Comments:
 *		This is the interface to all functions related to the solving stage of the computation:
 *			- compute_RLHS;
 *			- update_solution;
 *			- update_space.
 */

/*bool evaluate_exit_condition (const struct S_solver_info*const solver_info)
{
	if (solver_info->steady) {
		if (solver_info->linear)
			return true;

		double const EXIT_RHS_RATIO = 1e10;
//		double const maxRHS = compute_maxRHS();
double maxRHS = 0.0, maxRHS0 = 0.0;

		if ((maxRHS0/maxRHS > EXIT_RHS_RATIO) || (maxRHS < 1e1*EPS))
			return true;
	} else {
		if (time == final_time)
			return true;
	}
	return false;
}*/

void compute_final_solution (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Compute the final solution.
	 */

	for (bool finished = false; !finished; ) {
EXIT_UNSUPPORTED;
finished = false;
		compute_RLHS(solver_info);
//		finalize_RLHS(solver_info);
//		update_solution(solver_info);
//		update_space(solver_info); // For adaptation (ToBeDeleted: remove this comment)

//		finished = evaluate_exit_condition(solver_info);
	}

// Move to postprocessing function: (ToBeDeleted)
//	update_gradients();
//	if (solver_info->output)
//		output_to_paraview("SolFinal_");
}

void compute_RLHS (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Compute the (R)ight and optionally the (L)eft (H)and (S)ide residual contributions.
	 *
	 *	Comments:
	 *		The RHS includes all terms of the discretization except for the time-stepping term (i.e. d/dt What = RHS).
	 *		The LHS represents the linearization of the RHS (i.e. LHS = d(RHS)/d(What)). The LHS is then used to update
	 *		the solution using Newton's method:
	 *			1) approximate: RHS(What_exact) = 0 = RHS(What) + LHS(What)*dWhat + O(dWhat^2);
	 *			2) solve:       LHS*dWhat = -RHS(What);
	 *			3) update:      What += dWhat
	 *			4) repeat:      if LHS is a function of What, repeat until convergence.
	 */

	switch (solver_info->method) {
	case METHOD_DG: {
		compute_GradW_DG(solver_info);
		compute_VOLUME_RLHS_DG(solver_info);
		compute_FACE_RLHS_DG(solver_info);
		break;
	} case METHOD_HDG:
		compute_VOLUME_RLHS_HDG(solver_info);
		compute_FACE_RLHS_HDG(solver_info);
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

struct S_LHS_info constructor_LHS_info (double*const LHS, const struct S_VOLUME*const V0,
                                        const struct S_VOLUME*const V1, const InsertMode addv)
{
	/*
	 *	Purpose:
	 *		Initialize struct holding LHS data information.
	 *
	 *	Comments:
	 *		V0 gives the global row index and V1 the column index.
	 */

	struct S_LHS_info LHS_info;

	LHS_info.LHS     = LHS;
	LHS_info.VOLUME[0] = V0;
	LHS_info.VOLUME[1] = V1;
	LHS_info.IndA[0] = V0->IndA;
	LHS_info.IndA[1] = V1->IndA;
	LHS_info.Nn[0]   = V0->NvnS;
	LHS_info.Nn[1]   = V1->NvnS;
	LHS_info.addv    = addv;

	correct_collocated_for_symmetry_local(&LHS_info);

	return LHS_info;
}

void fill_PetscMat (const struct S_solver_info*const solver_info, const struct S_LHS_info*const LHS_info)
{
	/*
	 *	Purpose:
	 *		Add individual LHS contributions to the global system matrix.
	 */

	Mat A = solver_info->A;

	const double*const LHS = LHS_info->LHS;

	const size_t*const IndA = LHS_info->IndA,
	            *const Nn   = LHS_info->Nn;

	PetscInt idxm[Nn[0]],
	         idxn[Nn[1]];
	const InsertMode addv = LHS_info->addv;

	const unsigned int Neq  = DB.Neq,
	                   Nvar = DB.Nvar;

	for (size_t eq = 0; eq < Neq; eq++) {
		const size_t Indm = IndA[0] + eq*Nn[0];
		for (size_t i = 0; i < Nn[0]; i++)
			idxm[i] = Indm+i;

		for (size_t var = 0; var < Nvar; var++) {
			const size_t Indn = IndA[1] + var*Nn[1];
			for (size_t i = 0; i < Nn[1]; i++)
				idxn[i] = Indn+i;

			const PetscScalar*const vv = &LHS[(eq*Nvar+var)*Nn[0]*Nn[1]];

			MatSetValues(A,Nn[0],idxm,Nn[1],idxn,vv,addv);
		}
	}
}


void finalize_RLHS (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Finalize the (R)ight and (L) (H)and (S)ide terms.
	 *
	 *	Comments:
	 *		The original implementation uses redundant storage of RHS and LHS terms separated based on whether they are
	 *		computed using VOLUME or FACE integrals. This function combines these contributions if applicable.
	 *		(ToBeModified)
	 */

	switch (solver_info->method) {
	case METHOD_DG: {
// Refactor such that the RHS and LHS terms are directly store in the same memory location. (ToBeDeleted)
		EXIT_UNSUPPORTED;
		break;
	} case METHOD_HDG:
// Change info to RLHS (ToBeDeleted)
		EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}
}

void update_solution (const struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Update the solution using the appropriate explicit/implicit updating procedure.
	 *
	 *	Comments:
	 *		For the implicit solver, the LHS is then used to update the solution using Newton's method:
	 *			Taylor-Expansion: RHS(What_converged) := 0 = RHS(What) + LHS(What)*dWhat + O(dWhat^2).
	 *
	 *			1) solve:   LHS*dWhat = -RHS(What);
	 *			2) update:  What += dWhat;
	 *			3) iterate: if the LHS is a function of What (i.e. the equation is nonlinear), repeat until convergence.
	 */

	if (solver_info->imex_type == 'E') {
	} else if (solver_info->imex_type == 'I') {
/*		Mat A   = NULL;
		Vec b   = NULL,
		    x   = NULL;
		KSP ksp = NULL;

// Advection and Poisson (ToBeDeleted)
		solver_implicit_linear_system(&A,&b,&x,&ksp,0,PrintEnabled);
		solver_implicit_update_What(x);

		KSPDestroy(&ksp);
		finalize_ksp(&A,&b,&x,2);*/
	}
}

static void compute_nnz          (struct S_solver_info*const solver_info);
static void create_petsc_structs (struct S_solver_info*const solver_info);

void initialize_petsc_structs (struct S_solver_info*const solver_info)
{
	if (solver_info->imex_type != 'I')
		EXIT_UNSUPPORTED;

// Delete redundant functions in finalize_LHS when working (ToBeDeleted).
	set_global_indices(solver_info);
	compute_nnz(solver_info);
	create_petsc_structs(solver_info);
}

static void assemble_Mat (Mat A)
{
	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
}

static void assemble_Vec (Vec a)
{
	VecAssemblyBegin(a);
	VecAssemblyEnd(a);
}

void assemble_petsc_structs (struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Assemble Petsc data storage structs for global system solve.
	 *
	 *	Comments:
	 *		Called after containers are filled before passing to KSP.
	 */

	assemble_Mat(solver_info->A);
	if (solver_info->create_RHS) {
		assemble_Vec(solver_info->b);
		assemble_Vec(solver_info->x);
	}
}

void destroy_petsc_structs (struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Destroy Petsc data storage structs for global system solve.
	 */


	MatDestroy(&solver_info->A); solver_info->A = NULL;
	if (solver_info->create_RHS) {
		VecDestroy(&solver_info->b); solver_info->b = NULL;
		VecDestroy(&solver_info->x); solver_info->x = NULL;
	}
}



void set_global_indices (struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Set the indices of the degrees of freedom of VOLUME/FACE unknowns in the global system matrix.
	 *
	 *	Comments:
	 *		The total number of globally coupled degrees of freedom is also computed and stored.
	 *
	 *	Notation:
	 *		IndA_G  : (G)lobal (Ind)ex of global system matrix (A)
	 *		nnz_[1] : (n)umber of (n)on-(z)ero
	 *		          [1] : (d)iagonal, (o)ff-diagonal
	 */

	unsigned int Nvar = DB.Nvar;

	unsigned int IndA_G = 0;
	switch (solver_info->method) {
	case METHOD_DG:
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			VOLUME->IndA = IndA_G;
			IndA_G += Nvar*(VOLUME->NvnS);
// Currently needed for linearization testing (ToBeDeleted)
VOLUME->nnz_d = Nvar*(VOLUME->NvnS);
		}
		break;
	case METHOD_HDG:
		EXIT_UNSUPPORTED;
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			FACE->IndA = IndA_G;
			IndA_G += Nvar*(FACE->NvnS);
		}
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	solver_info->dof = IndA_G;
}

static void add_to_nnz (const unsigned int IndA, const unsigned int nnz_c, const unsigned int nnz_n,
                        PetscInt*const nnz);
static void add_to_nnz_FACE (const struct S_FACE*const FACE, const struct S_VOLUME*const VOLUME, PetscInt*const nnz);
static void update_VOLUME_FACE_pointers (void);

static void compute_nnz (struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Compute the (n)umber of (n)on-(z)ero entries in the global system matrix in each row.
	 */

	unsigned int Nvar = DB.Nvar;

	PetscInt* nnz = calloc((solver_info->dof) , sizeof *nnz); // keep

	switch (solver_info->method) {
	case METHOD_DG:
		// Diagonal contribution
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next)
			add_to_nnz(VOLUME->IndA,Nvar*(VOLUME->NvnS),Nvar*(VOLUME->NvnS),nnz);

		// Off-diagonal contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			struct S_VOLUME *VL = FACE->VL,
			                *VR = FACE->VR;

			if (VL != VR) {
				add_to_nnz(VL->IndA,Nvar*(VL->NvnS),Nvar*(VR->NvnS),nnz);
				add_to_nnz(VR->IndA,Nvar*(VR->NvnS),Nvar*(VL->NvnS),nnz);
			}
		}
		break;
	case METHOD_HDG:
		update_VOLUME_FACE_pointers();
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			// Diagonal contribution
			add_to_nnz(FACE->IndA,Nvar*(FACE->NvnS),Nvar*(FACE->NvnS),nnz);

			// Off-diagonal contributions
			struct S_VOLUME *VL = FACE->VL,
			                *VR = FACE->VR;

			add_to_nnz_FACE(FACE,VL,nnz);
			if (VL != VR)
				add_to_nnz_FACE(FACE,VR,nnz);
		}
		EXIT_UNSUPPORTED;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	solver_info->nnz = nnz;
}

static void create_petsc_structs (struct S_solver_info*const solver_info)
{
	/*
	 *	Purpose:
	 *		Initialize Petsc data storage structs for global system solve.
	 *
	 *	Comments:
	 *		Memory for the RHS and solution vector need not be allocated for the linearization tests.
	 */

	const unsigned int dof = solver_info->dof;
	PetscInt* nnz = solver_info->nnz;

	MatCreateSeqAIJ(MPI_COMM_WORLD,dof,dof,0,nnz,&solver_info->A); // keep
	if (solver_info->create_RHS) {
		VecCreateSeq(MPI_COMM_WORLD,dof,&solver_info->b); // keep
		VecCreateSeq(MPI_COMM_WORLD,dof,&solver_info->x); // keep
	}

	free(nnz); nnz = NULL;
}

static void add_to_nnz (const unsigned int IndA, const unsigned int nnz_c, const unsigned int nnz_n, PetscInt*const nnz)
{
	/*
	 *	Purpose:
	 *		Increase nnz based on the number of degrees of freedom in the (c)urrent and (n)eighbouring elements.
	 */

	for (size_t i = IndA; i < IndA+nnz_c; i++)
		nnz[i] += nnz_n;
}

static void add_to_nnz_FACE (const struct S_FACE*const FACE, const struct S_VOLUME*const VOLUME, PetscInt*const nnz)
{
	/*
	 *	Purpose:
	 *		Analogue of add_to_nnz where coupled degrees of freedom are between FACEs adjacent to neighbouring VOLUMEs.
	 */

	unsigned int Nvar = DB.Nvar;

	for (size_t f = 0; f < NFMAX; f++) {
	for (size_t lsf = 0; lsf < NSUBFMAX; lsf++) {
		struct S_FACE *FNeigh = VOLUME->FACE[f*NFMAX+lsf];
		if (!FNeigh)
			break;

		if (FACE == FNeigh)
			continue;

		add_to_nnz(FNeigh->IndA,Nvar*(FACE->NvnS),Nvar*(FNeigh->NvnS),nnz);
	}}
}

static void update_VOLUME_FACE_pointers (void)
{
	/*
	 *	Purpose:
	 *		Update VOLUME->FACE (the array of pointers to neighbouring FACEs for each VOLUME).
	 */

	for (struct S_VOLUME* VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		for (size_t i = 0; i < NFMAX*NSUBFMAX; i++)
			VOLUME->FACE[i] = NULL;
	}

	for (struct S_FACE* FACE = DB.FACE; FACE; FACE = FACE->next) {
		size_t const Nsides = (FACE->Boundary ? 1 : 2 );
		for (size_t i = 0; i < Nsides; i++)
			FACE->data[i].VOLUME->FACE[FACE->data[i].Indsfh] = FACE;
	}
}


// "class" functions

struct S_solver_info constructor_solver_info (const bool display, const bool output, const bool adapt,
                                              const char imex_type, const unsigned int method)
{
	struct S_solver_info solver_info;

	switch (DB.PDE_index) {
	case PDE_ADVECTION:
		solver_info.positivity = false;
		solver_info.symmetric  = false;
		solver_info.linear     = true;
		break;
	case PDE_POISSON:
		solver_info.positivity = false;
		solver_info.symmetric  = true;
		solver_info.linear     = true;
		break;
	case PDE_EULER:
	case PDE_NAVIERSTOKES:
		solver_info.positivity = true;
		solver_info.symmetric  = false;
		solver_info.linear     = false;
		break;
	default:
		EXIT_UNSUPPORTED;
		break;
	}

	solver_info.display = display;
	solver_info.output  = output;
	solver_info.adapt   = adapt;

	if (method != METHOD_DG && method != METHOD_HDG)
		EXIT_UNSUPPORTED;
	solver_info.method = method;

	if (imex_type != 'E' && imex_type != 'I')
		EXIT_UNSUPPORTED;
	solver_info.imex_type = imex_type;

	solver_info.compute_FACE = true;
	if (imex_type == 'I') {
		solver_info.create_RHS = true;

		solver_info.A = NULL;
		solver_info.b = NULL;
		solver_info.x = NULL;
	}

	return solver_info;
}
