// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Compute L2 errors for supported test cases.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnS, NvnI;
	double       *I_vG_vI, *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	P = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (1 || type == TRI || type == TET || type == PYR)
		ELEMENT_OPS = ELEMENT;
	else if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];

	OPS->NvnS = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGs_vIs[P][P][0];
		OPS->w_vI    = ELEMENT_OPS->w_vIs[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGc_vIc[P][P][0];
		OPS->w_vI    = ELEMENT_OPS->w_vIc[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];
	}
}

static void compute_errors_PeriodicVortex(void)
{
	// Initialize DB Parameters
	unsigned int d    = DB.d,
	             Nvar = DB.Nvar;
	double       MInf           = DB.MInf,
	             pInf           = DB.pInf,
	             TInf           = DB.TInf,
	             VInf           = DB.VInf,
	             Rg             = DB.Rg,
	             Cscale         = DB.Cscale,
	             PeriodL        = DB.PeriodL,
	             PeriodFraction = DB.PeriodFraction,
	             Xc0            = DB.Xc,
	             Yc             = DB.Yc,
	             Rc             = DB.Rc;

	// Standard datatypes
	unsigned int i, j, NvnS, NvnI, IndU;
	double       DistTraveled;
	double       Xc, uInf, vInf, rhoInf,
	             *XYZ_vI, *X_vI, *Y_vI, *r2, C,
	             *rho, *p, *s, *rhoEx, *uEx, *vEx, *wEx, *pEx, *sEx, *U, *UEx, *What, *W,
	             *detJV_vI, *w_vI, *ChiS_vI, *wdetJV_vI,
	             Vol, *L2Error, *L2Error2, err;

	struct S_OPERATORS *OPS;
	struct S_VOLUME *VOLUME;


	OPS = malloc(sizeof *OPS); // free

	L2Error  = malloc(NVAR3D+1 * sizeof *L2Error);  // free
	L2Error2 = calloc(NVAR3D+1 , sizeof *L2Error2); // free

	uInf    = MInf*sqrt(GAMMA*Rg*TInf);
	vInf   = 1.0;
	rhoInf = pInf/(Rg*TInf);
	C      = Cscale*VInf;

	DistTraveled = PeriodL*PeriodFraction;
	Xc = Xc0 + DistTraveled;
	while (Xc > 0.5*PeriodL)
		Xc -= PeriodL;

	Vol = 0;
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		init_ops(OPS,VOLUME,0);

		NvnS    = OPS->NvnS;
		NvnI    = OPS->NvnI;
		w_vI    = OPS->w_vI;
		ChiS_vI = OPS->ChiS_vI;

		detJV_vI = VOLUME->detJV_vI;

		wdetJV_vI = malloc(NvnI * sizeof *wdetJV_vI); // free
		for (i = 0; i < NvnI; i++)
			wdetJV_vI[i] = w_vI[i]*detJV_vI[i];

		W = malloc(NvnI*Nvar   * sizeof *W); // free
		U = malloc(NvnI*NVAR3D * sizeof *U); // free
		s = malloc(NvnI        * sizeof *s); // free

		What = VOLUME->What;
		mm_CTN_d(NvnI,Nvar,NvnS,ChiS_vI,What,W);

		convert_variables(W,U,d,3,NvnI,1,'c','p');

		rho = &U[NvnI*0];
		p   = &U[NvnI*(NVAR3D-1)];
		for (i = 0; i < NvnI; i++)
			s[i] = p[i]/pow(rho[i],GAMMA);

		UEx = malloc(NvnI*NVAR3D * sizeof *UEx); // free
		sEx = malloc(NvnI        * sizeof *sEx); // free

		rhoEx = &UEx[NvnI*0];
		uEx   = &UEx[NvnI*1];
		vEx   = &UEx[NvnI*2];
		wEx   = &UEx[NvnI*3];
		pEx   = &UEx[NvnI*4];

		XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
		mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

		X_vI = &XYZ_vI[0*NvnI];
		Y_vI = &XYZ_vI[1*NvnI];

		r2 = malloc(NvnI * sizeof *r2); // free
		for (i = 0; i < NvnI; i++)
			r2[i] = (pow(X_vI[i]-Xc,2.0)+pow(Y_vI[i]-Yc,2.0))/(Rc*Rc);

		for (i = 0; i < NvnI; i++) {
			uEx[i]   = uInf - C*(Y_vI[i]-Yc)/(Rc*Rc)*exp(-0.5*r2[i]);
			vEx[i]   = vInf + C*(X_vI[i]-Xc)/(Rc*Rc)*exp(-0.5*r2[i]);
			wEx[i]   = 0.0;
			pEx[i]   = pInf - rhoInf*(C*C)/(2*Rc*Rc)*exp(-r2[i]);
			rhoEx[i] = rhoInf;
			sEx[i]   = pEx[i]/pow(rhoEx[i],GAMMA);
		}

		for (i = 0; i <= NVAR3D; i++) {
			IndU = i*NvnI;
			if (i == 0 || i == NVAR3D-1) { // rho, p
				for (j = 0; j < NvnI; j++) {
					err = (U[IndU+j]-UEx[IndU+j])/UEx[IndU+j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else if (i > 0 && i < NVAR3D-1) { // u, v, w (Not normalized as variables may be negative)
				for (j = 0; j < NvnI; j++) {
					err = U[IndU+j]-UEx[IndU+j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else if (i == NVAR3D) { // s
				for (j = 0; j < NvnI; j++) {
					err = (s[j]-sEx[j])/sEx[j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			}
		}
	
		for (i = 0; i < NvnI; i++)
			Vol += wdetJV_vI[i];

		free(wdetJV_vI);
		free(W);
		free(U);
		free(s);
		free(UEx);
		free(sEx);
		free(XYZ_vI);
		free(r2);
	}
	free(OPS);

	for (i = 0; i <= NVAR3D; i++)
		L2Error[i] = sqrt(L2Error2[i]/Vol);

	free(L2Error2);

	// Write to file
array_print_d(1,NVAR3D+1,L2Error,'R');


	free(L2Error);
}

void compute_errors(void)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	if (strstr(TestCase,"dSphericalBump") != NULL)
		; // compute_errors_dSphericalBump();
	else if (strstr(TestCase,"GaussianBump") != NULL)
		; // compute_errors_GaussianBump();
	else if (strstr(TestCase,"PeriodicVortex") != NULL)
		compute_errors_PeriodicVortex();
	else if (strstr(TestCase,"PolynomialBump") != NULL)
		; // compute_errors_PolynomialBump();
	else if (strstr(TestCase,"SupersonicVortex") != NULL)
		; // compute_errors_SupersonicVortex();
}
