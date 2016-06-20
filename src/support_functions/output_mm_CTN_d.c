// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*
 *	Purpose:
 *		Output fully unrolled code for mv multiplications.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

#define STRLEN_MIN 60
#define STRLEN_MAX 508

static void fprintf_t(FILE *fID, unsigned int Ntabs, unsigned int includeNewLine, const char *String)
{
	unsigned int i;

	for (i = 0; i < Ntabs; i++)
		fprintf(fID,"\t");
	fprintf(fID,"%s",String);
	if (includeNewLine)
		fprintf(fID,"\n");
}

int main(void)
{
	char         f_name[STRLEN_MAX];
	unsigned int i, j, m, k, mMax_unrolled, kMax_unrolled;

	FILE *fID;

	strcpy(f_name,"mm_CTN_d_unrolled_mv.txt");

	mMax_unrolled = 12;
	kMax_unrolled = 12;


	if ((fID = fopen(f_name,"w")) == NULL)
		printf("Error: File: %s, did not open.\n",f_name), exit(1);

	fprintf_t(fID,2,1,"switch(m) {");
	for (m = 1; m <= mMax_unrolled; m++) {
		if (m == 1)
			fprintf(fID,"\t\tcase %d: {\n",m);
		else
			fprintf(fID,"\t\t} case %d: {\n",m);

		fprintf_t(fID,3,1,"switch(k) {");
		for (k = 1; k <= kMax_unrolled; k++) {
			if (k == 1)
				fprintf(fID,"\t\t\tcase %d: {\n",k);
			else
				fprintf(fID,"\t\t\t} case %d: {\n",k);

			// Initialize pointers to A
			for (i = 0; i < m; i++) {
				if (!i)
					fprintf_t(fID,4,0,"register double");
				else
					fprintf_t(fID,4,0,"               ");

				for (j = 0; j < k; j++) {
					if (!i && !j)
						fprintf(fID," *a%-2d = A   ",i*k+j);
					else
						fprintf(fID," *a%-2d = A+%-2d",i*k+j,i*k+j);
					fprintf(fID,",");
				}
				fprintf(fID,"\n");
			}

			// Initialize pointers to B
			fprintf_t(fID,4,0,"               ");
			for (j = 0; j < k; j++) {
				if (!j)
					fprintf(fID," *b%-2d = B   ",j);
				else
					fprintf(fID," *b%-2d = B+%-2d",j,j);

				if (j < k-1 && k != 1)
					fprintf(fID,",");
				else
					fprintf(fID,";");
			}
			fprintf(fID,"\n\n");

			// Print unrolled mv computation of C = A*B

			// First column
			for (i = 0; i < m; i++) {
				if (!i)
					fprintf_t(fID,4,0,"*C    ");
				else
					fprintf_t(fID,4,0,"*(++C)");
				for (j = 0; j < k; j++) {
					if (!j)
						fprintf(fID," = (*a%-2d)*(*b%-2d)",i*k+j,j);
					else
						fprintf(fID," + (*a%-2d)*(*b%-2d)",i*k+j,j);
				}
				fprintf(fID,";\n");
			}
			// Other columns
			fprintf(fID,"\n");
			fprintf_t(fID,4,1,"for (register unsigned int nRem = n-1; nRem--; ) {");
			fprintf_t(fID,5,0,"");
			for (j = 0; j < k; j++) {
				fprintf(fID,"b%-2d += k",j);
				if (j != k-1)
					fprintf(fID,", ");
			}
			fprintf(fID,";\n\n");
			for (i = 0; i < m; i++) {
				fprintf_t(fID,5,0,"*(++C)");
				for (j = 0; j < k; j++) {
					if (!j)
						fprintf(fID," = (*a%-2d)*(*b%-2d)",i*k+j,j);
					else
						fprintf(fID," + (*a%-2d)*(*b%-2d)",i*k+j,j);
				}
				fprintf(fID,";\n");
			}
			fprintf_t(fID,4,1,"}");
			fprintf(fID,"\t\t\t\tbreak; // m%dk%d\n",m,k);
		}
		fprintf_t(fID,3,1,"} default: {");
		fprintf(fID,"\t\t\t\tprintf(\"Error: Unsupported m = %%d, k = %%d, combination in mm_CTN_d (useBlas = 0).\\n\",m,k), exit(1);\n");
		fprintf(fID,"\t\t\t\tbreak; // m%d\n",m);
		fprintf_t(fID,3,1,"}}");

		fprintf_t(fID,3,1,"break;");
	}
	fprintf_t(fID,2,1,"} default: {");
	fprintf(fID,"\t\t\tprintf(\"Error: Unsupported m = %%d, in mm_CTN_d (useBlas = 0).\\n\",m), exit(1);\n");
	fprintf(fID,"\t\t\tbreak;\n");
	fprintf_t(fID,2,1,"}}");

	fclose(fID);

	return 0;
}
