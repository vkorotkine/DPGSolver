// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/*
 *	Purpose:
 *		Output convergence orders of the selected test case in tabular format.
 *
 *	Comments:
 *		table_to_latex.c can then be used to convert output from this function to a latex compatible format.
 *
 *	Notation:
 *
 *	References:
 *
 */

#define STRLEN_MIN 60
#define STRLEN_MAX 508

#define EPS        1.0e-15
#define DMAX       3

static void table_to_latex(const unsigned int d, const unsigned int NVars, const unsigned int MLMin,
                           const unsigned int MLMax, const unsigned int PMin, const unsigned int PMax,
                           const unsigned int *CasesRun, const double *h, double **L2Errors, double **ConvOrders,
                           char *TestCase, char *MeshType);

int main(void)
{
	unsigned int Testing;

	char         *TestCase, *MeshType, f_name[STRLEN_MAX], string[STRLEN_MIN], *data, StringRead[STRLEN_MAX];
	unsigned int i, ML, P,
	             Indh,
	             d, NVars, MLMin, MLMax, PMin, PMax, NML, NP;
	int          offset;
	double       **L2Errors, **ConvOrders, *h, tmp_d;

	FILE *fID;

	TestCase = malloc(STRLEN_MAX * sizeof *TestCase); // free
	MeshType = malloc(STRLEN_MAX * sizeof *MeshType); // free

	Testing = 1;

	strcpy(TestCase,"PeriodicVortex");
	strcpy(TestCase,"SupersonicVortex");
	strcpy(TestCase,"Poisson");
//	strcpy(MeshType,"StructuredTRI");
	strcpy(MeshType,"ToBeCurvedStructuredTRI");
//	strcpy(MeshType,"ToBeCurvedStructuredQUAD");
//	strcpy(MeshType,"ToBeCurvedStructuredTET");
//	strcpy(MeshType,"ToBeCurvedStructuredHEX");

	d     = 2;
//	NVars = DMAX+2+1;
	NVars = DMAX+1;
	MLMin = 1; MLMax = 5; NML = MLMax-MLMin+1;
	PMin  = 1; PMax  = 3; NP  = PMax-PMin+1;

	unsigned int CasesRun[45] = { 1, 1, 1, 1, 1,
	                              1, 1, 1, 1, 1,
	                              1, 1, 1, 1, 1,
	                              1, 1, 1, 1, 1,
	                              1, 1, 1, 1, 1,
	                              1, 1, 1, 1, 1,
	                              1, 1, 0, 0, 0,
	                              0, 0, 0, 0, 0,
	                              0, 0, 0, 0, 0};

	if (Testing) {
		for (ML = 0; ML <= MLMax; ML++) {
		for (P = 0;  P  <= PMax;  P++) {
			Indh = ML*NP+P;
			CasesRun[Indh] = 1;
		}}
	}

	L2Errors   = malloc(NVars * sizeof *L2Errors);   // free
	ConvOrders = malloc(NVars * sizeof *ConvOrders); // free
	for (i = 0; i < NVars; i++) {
		L2Errors[i]   = calloc(NML*NP , sizeof *L2Errors[i]);   // free
		ConvOrders[i] = calloc(NML*NP , sizeof *ConvOrders[i]); // free
	}
	h = calloc(NML*NP , sizeof *h); // free

	// Read in data and compute convergence orders
	for (ML = MLMin; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);
		if (CasesRun[Indh]) {
			strcpy(f_name,"../../cases/results/");
			strcat(f_name,TestCase); strcat(f_name,"/");
			strcat(f_name,MeshType); strcat(f_name,"/");
			strcat(f_name,"L2errors_");
			sprintf(string,"%dD_",d);   strcat(f_name,string);
										strcat(f_name,MeshType);
			sprintf(string,"_ML%d",ML); strcat(f_name,string);
			sprintf(string,"P%d",P);    strcat(f_name,string);
			strcat(f_name,".txt");

			if ((fID = fopen(f_name,"r")) == NULL)
				printf("Error: File: %s, did not open.\n",f_name), exit(1);

			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				i = 0;
				data = StringRead;
				if (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
					data += offset;
					h[Indh] = 1.0/pow(tmp_d,1.0/d);
				}
				while (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
					L2Errors[i++][Indh] = tmp_d;
					data += offset;
				}
			}
			fclose(fID);
		}
	}}

	for (ML = MLMin+1; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);
		if (CasesRun[Indh]) {
			for (i = 0; i < NVars; i++) {
				if (fabs(h[Indh]) > EPS && fabs(h[Indh-NP]) > EPS)
					ConvOrders[i][Indh] = log10(L2Errors[i][Indh]/L2Errors[i][Indh-NP])/log10(h[Indh]/h[Indh-NP]);
			}
		}
	}}

	// Output to file
	table_to_latex(d,NVars,MLMin,MLMax,PMin,PMax,CasesRun,h,L2Errors,ConvOrders,TestCase,MeshType);

	free(TestCase);
	free(MeshType);

	for (i = 0; i < NVars; i++) {
		free(L2Errors[i]);
		free(ConvOrders[i]);
	}
	free(L2Errors);
	free(ConvOrders);
	free(h);

	return 0;
}

static void table_to_latex(const unsigned int d, const unsigned int NVars, const unsigned int MLMin,
                           const unsigned int MLMax, const unsigned int PMin, const unsigned int PMax,
                           const unsigned int *CasesRun, const double *h, double **L2Errors, double **ConvOrders,
                           char *TestCase, char *MeshType)
{
	char         **Vars_c, caption[STRLEN_MAX];
	unsigned int i, j, ML, P, NP, NVarsOut, Indp, Indh, IndVars[6], P_Print;

	FILE *fID;

	if (strstr(TestCase,"PeriodicVortex") != NULL ||
		strstr(TestCase,"SupersonicVortex") != NULL) {
		if      (d == 2) Indp = 3;
		else if (d == 3) Indp = 4;

		NVarsOut = NVars+d-3;
		Vars_c = malloc(NVarsOut * sizeof *Vars_c); // free
		for (i = 0; i < NVarsOut; i++)
			Vars_c[i] = malloc(STRLEN_MIN * sizeof *Vars_c[i]); // free

		strcpy(Vars_c[0],"$\\rho$");
		strcpy(Vars_c[1],"$u$    ");
		strcpy(Vars_c[2],"$v$    ");
		strcpy(Vars_c[Indp],"$p$    ");
		strcpy(Vars_c[Indp+1],"$s$    ");
		if (d == 3)
			strcpy(Vars_c[3],"$w$    ");

		for (i = 0; i < NVarsOut; i++) {
			if (d == 3 || i < Indp)
				IndVars[i] = i;
			else
				IndVars[i] = i+1;
		}
	} else if (strstr(TestCase,"Poisson")) {
		NVarsOut = NVars+d-3;
		Vars_c = malloc(NVarsOut * sizeof *Vars_c); // free
		for (i = 0; i < NVarsOut; i++)
			Vars_c[i] = malloc(STRLEN_MIN * sizeof *Vars_c[i]); // free

		strcpy(Vars_c[0],"$u$    ");
		strcpy(Vars_c[1],"$q_1$  ");
		strcpy(Vars_c[2],"$q_2$  ");
		if (d == 3)
			strcpy(Vars_c[3],"$q_3$  ");

		for (i = 0; i < NVarsOut; i++)
			IndVars[i] = i;
	}
	strcpy(caption,"Errors and Convergence Orders - ");
	strcat(caption,MeshType);
	strcat(caption," meshes");

	if ((fID = fopen("L2errs+Convergence.txt","w")) == NULL)
		printf("Error: File in table_to_latex did not open.\n"), exit(1);

	fprintf(fID,"\\begin{table}[!htbp]\n");
	fprintf(fID,"\\begin{center}\n");
	fprintf(fID,"\\caption{ %s }\n",caption);
	fprintf(fID,"\\resizebox{\\textwidth}{!}{\n");
	fprintf(fID,"\\begin{tabular}{| l | l | ");
	for (i = 0; i < 2; i++) {
		for (j = 0; j < NVarsOut; j++)
			fprintf(fID,"c ");
		fprintf(fID,"| ");
	}
	fprintf(fID,"}\n");

	fprintf(fID,"\t\\hline\n");
	fprintf(fID,"\t & & ");
	for (i = 0; i < 2; i++) {
		if      (i == 0) fprintf(fID," $L_2$ Error ");
		else if (i == 1) fprintf(fID," Conv. Order ");
		for (j = 0; j < NVarsOut; j++) {
			if (!(i == 1 && j == NVarsOut-1))
				fprintf(fID,"& ");
		}
	}
	fprintf(fID,"\\\\\n");
	fprintf(fID,"\t\\hline\n");
	fprintf(fID,"\tP. Order & Mesh Size ");
	for (i = 0; i < 2; i++) {
	for (j = 0; j < NVarsOut; j++) {
		fprintf(fID,"& %s ",Vars_c[j]);
	}}
	fprintf(fID,"\\\\\n");

	NP = PMax-PMin+1;
	for (P = PMin; P <= PMax; P++) {
		P_Print = 1;
		fprintf(fID,"\t\\hline\n");
		for (ML = MLMin; ML <= MLMax; ML++) {
			Indh = (ML-MLMin)*NP+(P-PMin);
			if (CasesRun[Indh]) {
				if (P_Print) {
					P_Print = 0;
					fprintf(fID,"P%1d\t& % .3e",P,h[Indh]);
					for (i = 0; i < 2; i++) {
					for (j = 0; j < NVarsOut; j++) {
						if      (i == 0) fprintf(fID," & % .3e",L2Errors[IndVars[j]][Indh]);
						else if (i == 1) fprintf(fID," & -");
					}}
					fprintf(fID," \\\\\n");
				} else {
					fprintf(fID,"\t& % .3e",h[Indh]);
					for (i = 0; i < 2; i++) {
					for (j = 0; j < NVarsOut; j++) {
						if      (i == 0) fprintf(fID," & % .3e",L2Errors[IndVars[j]][Indh]);
						else if (i == 1) fprintf(fID," & % .3e",ConvOrders[IndVars[j]][Indh]);
					}}
					fprintf(fID," \\\\\n");
				}
			}
		}
	}

	fprintf(fID,"\t\\hline\n");
	fprintf(fID,"\\end{tabular}\n");
	fprintf(fID,"}\n");
	fprintf(fID,"\\end{center}\n");
	fprintf(fID,"\\end{table}");

	for (i = 0; i < NVarsOut; i++)
		free(Vars_c[i]);
	free(Vars_c);

	fclose(fID);
}
