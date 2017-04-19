// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "explicit_VOLUME_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"

#include "solver_functions.h"
#include "sum_factorization.h"
#include "matrix_functions.h" // ToBeDeleted
#include "fluxes_inviscid.h"
#include "fluxes_viscous.h"
#include "array_free.h"
#include "array_print.h"

/*
 *	Purpose:
 *		Evaluate the VOLUME contributions to the RHS term.
 *
 *	Comments:
 *		Vectorization does not improve performance based on preliminary testing. This is likely a result of the large
 *		memory allocation/deallocation overhead which is required for the adaptive code. For non-adaptive versions of
 *		the code, this flexibility is not required and comparison may yield favourable results for the vectorized code.
 *		However, the ultimate goal of the code is only to run in the adaptive setting, thus this may not be worth
 *		pursuing. (ToBeModified)
 *
 *		EFE stands for (E)xact (F)lux (E)valuation, meaning that the flux is not represented as a polynomial and then
 *		interpolated to the cubature nodes. EFE is the analogue of the (C)hain(R)ule approach for the strong form of the
 *		scheme. See Zwanenburg(2016) for additional discussion.
 *
 *	Notation:
 *
 *	References:
 *		Zwanenburg(2016)-Equivalence_between_the_Energy_Stable_Flux_Reconstruction_and_Discontinuous_Galerkin_Schemes
 */

static void compute_Inviscid_VOLUME_RHS_EFE (void);
static void compute_Viscous_VOLUME_RHS_EFE  (void);
static void compute_VOLUMEVec_RHS_EFE       (void);

void explicit_VOLUME_info(void)
{
	if (DB.EFE) {
		switch (DB.Vectorized) {
		case 0:
			compute_Inviscid_VOLUME_RHS_EFE();
			compute_Viscous_VOLUME_RHS_EFE();
			break;
		default:
			compute_VOLUMEVec_RHS_EFE();
			break;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

static void compute_Inviscid_VOLUME_RHS_EFE(void)
{
	// Initialize DB Parameters
	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   Neq  = d+2;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI
			unsigned int const NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI = VOLUME->What;
			} else {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // free
				coef_to_values_vI(VDATA,'W');
			}

			// Compute Flux in reference space
			double *const F_vI = malloc(NvnI*d*Neq * sizeof *F_vI); // free

			flux_inviscid(NvnI,1,VDATA->W_vI,F_vI,d,Neq);

			if (!DB.Collocated)
				free(VDATA->W_vI);

			// Convert to reference space
			double *const Fr_vI = malloc(NvnI*Neq*d * sizeof *Fr_vI); // free
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			// Compute RHS term
			unsigned int const NvnS = VDATA->OPS[0]->NvnS;

			memset(VOLUME->RHS,0.0,NvnS*Neq * sizeof *(VOLUME->RHS));
			finalize_VOLUME_Inviscid_Weak(Neq,Fr_vI,VOLUME->RHS,'E',VDATA);
			free(Fr_vI);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}

static void compute_Viscous_VOLUME_RHS_EFE(void)
{
	/*
	 *	Purpose:
	 *		Add contributions from the viscous term to the RHS.
	 *
	 *	Comments:
	 *		The viscous VOLUME contributions have a nearly identical form to those of the inviscid contributions.
	 *		Consider combining the two functions in the future. (ToBeModified)
	 */

	if (!DB.Viscous)
		return;

	unsigned int const d    = DB.d,
	                   Nvar = d+2,
	                   Neq  = d+2;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	if (strstr(DB.Form,"Weak")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			init_VDATA(VDATA,VOLUME);

			// Obtain W_vI and Q_vI
			unsigned int const NvnI = VDATA->OPS[0]->NvnI;
			if (DB.Collocated) {
				VDATA->W_vI = VOLUME->What;
				VDATA->Q_vI = VOLUME->Qhat;
			} else {
				VDATA->W_vI = malloc(NvnI*Nvar * sizeof *(VDATA->W_vI)); // free
				VDATA->Q_vI = malloc(d         * sizeof *(VDATA->Q_vI)); // free
				for (size_t dim = 0; dim < d; dim++)
					VDATA->Q_vI[dim] = malloc(NvnI*Nvar * sizeof *(VDATA->Q_vI[dim])); // free

				coef_to_values_vI(VDATA,'W');
				coef_to_values_vI(VDATA,'Q');
			}

			// Compute negated Flux in reference space
			double *const F_vI = malloc(NvnI*d*Neq * sizeof *F_vI);

			flux_viscous(NvnI,1,VDATA->W_vI,(const double *const *const) VDATA->Q_vI,F_vI);

			if (!DB.Collocated) {
				free(VDATA->W_vI);
				array_free2_d(d,VDATA->Q_vI);
			}

			// Convert to reference space
			double *const Fr_vI = malloc(NvnI*Neq*d * sizeof *Fr_vI); // free
			convert_between_rp(NvnI,Neq,VOLUME->C_vI,F_vI,Fr_vI,"FluxToRef");
			free(F_vI);

			// Compute RHS term
			finalize_VOLUME_Viscous_Weak(Neq,Fr_vI,VOLUME->RHS,'E',VDATA);
			free(Fr_vI);
		}
	} else if (strstr(DB.Form,"Strong")) {
		EXIT_UNSUPPORTED;
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
}


#undef I // No complex variables used here
static void compute_VOLUMEVec_RHS_EFE(void)
{
	// POTENTIALLY INEFFICIENT DUE TO LACK OF CODING EXPERIENCE WHEN THIS FUNCTION WAS WRITTEN. REVISIT. ToBeDeleted

	// Initialize DB Parameters
	char         *Form = DB.Form;
	unsigned int d          = DB.d,
	             Collocated = DB.Collocated,
	             Nvar       = DB.Nvar,
	             Neq        = DB.Neq,
	             NP         = DB.NP,
	             NECgrp     = DB.NECgrp,
	             *NVgrp     = DB.NVgrp,
	             ***SF_BE   = DB.SF_BE;

	// Standard datatypes
	unsigned int i, j, k, eq, dim1, dim2, P, curved, Eclass, iMax, jMax,
	             IndVgrp, IndC, IndF, IndFr, IndRHS, NvnI, NvnS, NV,
	             NvnS_SF[2], NvnI_SF[2], NIn[3], NOut[3], Diag[3];
	double       *What_vS, *W_vI, **What_vS_ptr, *What, *WhatVec, *F_vI, *Fr_vI,
	             **C_vI_ptr, *C_vIl, *C_vI, *C_vIVec,
	             **RHS_ptr, *RHSl, *RHS, *RHSVec, *DFr;
	double const *OP0, *OP1, *OP[3];

	struct S_OpCSR const *const *D_sp;

	struct S_OPERATORS_V *OPS[2];
	struct S_VOLUME      *VOLUME;

	for (i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	What_vS_ptr = malloc(Nvar * sizeof *What_vS_ptr); // free
	C_vI_ptr    = malloc(d*d  * sizeof *C_vI_ptr);    // free
	RHS_ptr     = malloc(Neq  * sizeof *RHS_ptr);     // free

	if (strstr(Form,"Weak")) {
		for (Eclass = 0; Eclass < NECgrp; Eclass++) {
		for (P = 0; P < NP; P++) {
		for (curved = 0; curved < 2; curved++) {
			IndVgrp = Eclass*NP*2 + P*2 + curved;
			VOLUME = DB.Vgrp[IndVgrp];

			if (VOLUME == NULL)
				continue;

			NV = NVgrp[IndVgrp];

			// Obtain operators
			init_ops_VOLUME(OPS[0],VOLUME,0);
			if (VOLUME->type == WEDGE)
				init_ops_VOLUME(OPS[1],VOLUME,1);

			// Obtain Vectorized W_vI
			NvnI = OPS[0]->NvnI;
			NvnS = OPS[0]->NvnS;
			W_vI = malloc(NvnI*NV*Nvar * sizeof *W_vI); // free
			C_vI = malloc(NvnI*NV*d*d  * sizeof *C_vI); // free

			if (Collocated) {
				for (i = 0; i < Nvar; i++)
					What_vS_ptr[i] = &W_vI[i*NV*NvnI];
			} else {
				What_vS = malloc(NvnS*NV*Nvar * sizeof *What_vS); // free

				for (i = 0; i < Nvar; i++)
					What_vS_ptr[i] = &What_vS[i*NV*NvnS];
			}
			for (i = 0, iMax = d*d; i < iMax; i++)
				C_vI_ptr[i] = &C_vI[i*NV*NvnI];

			for (i = 0, iMax = NV; i < iMax; i++) {
				if (i)
					VOLUME = VOLUME->grpnext;

				What = VOLUME->What;
				for (j = 0; j < Nvar; j++) {
					WhatVec = What_vS_ptr[j];
					for (k = 0; k < NvnS; k++)
						*WhatVec++ = *What++;
					What_vS_ptr[j] = WhatVec;
				}

				C_vIl = VOLUME->C_vI;
				for (j = 0, jMax = d*d; j < jMax; j++) {
					C_vIVec = C_vI_ptr[j];
					for (k = 0; k < NvnI; k++)
						*C_vIVec++ = *C_vIl++;
					C_vI_ptr[j] = C_vIVec;
				}
			}

			if (!Collocated) {
				if (Eclass == C_TP && SF_BE[P][0][0]) {
					for (i = 0; i < 1; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,0,0,NULL,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(What_vS,W_vI,NIn,NOut,NV*Nvar,OP,Diag,d);
				} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
					for (i = 0; i < 2; i++) {
						NvnS_SF[i] = OPS[i]->NvnS_SF;
						NvnI_SF[i] = OPS[i]->NvnI_SF;
					}
					get_sf_parameters(NvnS_SF[0],NvnI_SF[0],OPS[0]->ChiS_vI,
					                  NvnS_SF[1],NvnI_SF[1],OPS[1]->ChiS_vI,NIn,NOut,OP,d,3,Eclass);

					for (dim2 = 0; dim2 < d; dim2++)
						Diag[dim2] = 0;

					sf_apply_d(What_vS,W_vI,NIn,NOut,NV*Nvar,OP,Diag,d);
				} else {
					mm_CTN_d(NvnI,NV*Nvar,NvnS,OPS[0]->ChiS_vI,What_vS,W_vI);
				}

				free(What_vS);
			}

//array_print_d(NvnI*NV,Nvar,W_vI,'C');

			// Compute Flux in reference space
			F_vI = malloc(NvnI*NV*d*Neq * sizeof *F_vI); // free
			flux_inviscid(NvnI,NV,W_vI,F_vI,d,Neq);
			free(W_vI);

//array_print_d(NvnI*NV,Nvar*d,F_vI,'C');

			Fr_vI = calloc(NvnI*NV*d*Neq , sizeof *Fr_vI); // free
			for (eq = 0; eq < Neq; eq++) {
			for (dim1 = 0; dim1 < d; dim1++) {
			for (dim2 = 0; dim2 < d; dim2++) {
				IndFr = (eq*d+dim1)*NvnI*NV;
				IndF  = (eq*d+dim2)*NvnI*NV;
				IndC  = (dim1*d+dim2)*NvnI*NV;
				for (i = 0, iMax = NvnI*NV; i < iMax; i++)
					Fr_vI[IndFr+i] += F_vI[IndF+i]*C_vI[IndC+i];
			}}}
			free(F_vI);
			free(C_vI);

//array_print_d(NvnI*NV,Nvar*d,Fr_vI,'C');

			// Compute RHS terms
			DFr = malloc(NvnS*NV * sizeof *DFr); // free
			RHS = calloc(NvnS*NV*Neq , sizeof *RHS); // free
			if (Eclass == C_TP && SF_BE[P][0][0]) {
				for (i = 0; i < 1; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				double const *const        I = OPS[0]->I_Weak;
				double const *const *const D = OPS[0]->D_Weak;

				for (dim1 = 0; dim1 < d; dim1++) {
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],I,NvnI_SF[0],NvnS_SF[0],D[0],NIn,NOut,OP,d,dim1,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr,NIn,NOut,NV,OP,Diag,d);

						IndRHS = eq*NvnS*NV;
						for (i = 0, iMax = NvnS*NV; i < iMax; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Eclass == C_WEDGE && SF_BE[P][1][0]) {
				for (i = 0; i < 2; i++) {
					NvnS_SF[i] = OPS[i]->NvnS_SF;
					NvnI_SF[i] = OPS[i]->NvnI_SF;
				}

				for (dim1 = 0; dim1 < d; dim1++) {
					if (dim1 < 2) OP0 = OPS[0]->D_Weak[dim1], OP1 = OPS[1]->I_Weak;
					else          OP0 = OPS[0]->I_Weak,       OP1 = OPS[1]->D_Weak[0];
					get_sf_parameters(NvnI_SF[0],NvnS_SF[0],OP0,NvnI_SF[1],NvnS_SF[1],OP1,NIn,NOut,OP,d,3,Eclass);

					if (Collocated) {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 2;
						if (dim1 < 2)
							Diag[0] = 0;
						else
							Diag[dim1] = 0;
					} else {
						for (dim2 = 0; dim2 < d; dim2++)
							Diag[dim2] = 0;
					}

					for (eq = 0; eq < Neq; eq++) {
						sf_apply_d(&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr,NIn,NOut,NV,OP,Diag,d);

						IndRHS = eq*NvnS*NV;
						for (i = 0, iMax = NvnS*NV; i < iMax; i++)
							RHS[IndRHS+i] += DFr[i];
					}
				}
			} else if (Collocated && (Eclass == C_TP || Eclass == C_WEDGE)) {
				D_sp = OPS[0]->D_Weak_sp;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_CSR_d(NvnS,NV,NvnI,D_sp[dim1],&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr);

					IndRHS = eq*NvnS*NV;
					for (i = 0, iMax = NvnS*NV; i < iMax; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			} else {
				DFr = malloc(NvnS*NV * sizeof *DFr); // free
				double const *const *const D = OPS[0]->D_Weak;

				for (eq = 0; eq < Neq; eq++) {
				for (dim1 = 0; dim1 < d; dim1++) {
					mm_CTN_d(NvnS,NV,NvnI,D[dim1],&Fr_vI[(eq*d+dim1)*NvnI*NV],DFr);

					IndRHS = eq*NvnS*NV;
					for (i = 0, iMax = NvnS*NV; i < iMax; i++)
						RHS[IndRHS+i] += DFr[i];
				}}
			}
			free(DFr);
			free(Fr_vI);
//array_print_d(NvnS*NV,Neq,RHS,'C');

			// Copy values of Global RHS into local VOLUME RHS arrays
			// Profile and, if expensive, keep global RHS array
			for (i = 0; i < Neq; i++)
				RHS_ptr[i] = &RHS[i*NV*NvnS];

			for (VOLUME = DB.Vgrp[IndVgrp]; VOLUME; VOLUME = VOLUME->grpnext) {
				VOLUME->RHS = malloc(NvnS*Neq * sizeof *RHSl); // keep
				RHSl = VOLUME->RHS;

				for (eq = 0; eq < Neq; eq++) {
					RHSVec = RHS_ptr[eq];
					for (i = 0; i < NvnS; i++)
						*RHSl++ = *RHSVec++;
					RHS_ptr[eq] = RHSVec;
				}
//array_print_d(NvnS,Neq,VOLUME->RHS,'C');
			}
			free(RHS);
//exit(1);

		}}}
	}

	free(What_vS_ptr);
	free(C_vI_ptr);
	free(RHS_ptr);

	for (i = 0; i < 2; i++)
		free(OPS[i]);
}
