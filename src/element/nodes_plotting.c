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
/// \file

#include "nodes_plotting.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_elements.h"
#include "definitions_nodes.h"
#include "definitions_tol.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

// Static function declarations ************************************************************************************* //

///\{ \name Values of VTK linear cell types. See Figure 2 in $ROOT/doc/VTK_file_formats.pdf.
#define VTK_LINE       3
#define VTK_TRIANGLE   5
#define VTK_QUAD       9
#define VTK_TETRA      10
#define VTK_HEXAHEDRON 12
#define VTK_WEDGE      13
#define VTK_PYR        14

#define VTK_POLY_LINE  4
///\}

/** \brief Sets up (and allocates memory for) coordinates and connectivity information of plotting nodes.
 *  This function needs modernizing.
 */
static void plotting_element_info
	(double **rst,    ///< Data of \ref Nodes::rst.
	 int **connect,   ///< Data of \ref Plotting_Nodes::connect.
	 int **types,     ///< Data of \ref Plotting_Nodes::vtk_types.
	 int **connect_e, ///< Data of \ref Plotting_Nodes::connect_e.
	 int *n_n,        ///< Number of nodes (rst->ext_0).
	 int *n_e,        ///< Number of sub-elements (connect->ext_0);
	 int *d_,         ///< \ref Element::d.
	 int *n_edge,     ///< \ref Element::n_e.
	 int *n_n_edge,   ///< Number of nodes on each edge.
	 const int p,     ///< Order of the plotting nodes.
	 const int e_type ///< \ref Element::type.
	);

/** \brief Get the number of corners of the element of the input vtk_type.
 *  \return See brief. */
static int get_vtk_n_corners
	(const int vtk_type ///< The vtk element type.
	);

/// \brief `mutable` version of \ref destructor_const_Plotting_Nodes.
static void destructor_Plotting_Nodes
	(struct Plotting_Nodes* p_nodes ///< Standard.
	);

/// \brief `mutable` version of \ref destructor_const_Plotting_Nodes_part.
static void destructor_Plotting_Nodes_part
	(struct Plotting_Nodes* p_nodes ///< Standard.
	);

// Constructor functions ******************************************************************************************** //

const struct const_Plotting_Nodes* constructor_const_Plotting_Nodes (const int p, const int e_type)
{
	double* rst    = NULL;
	int* connect   = NULL,
	   * connect_e = NULL,
	   * vtk_types = NULL;
	int n_n      = -1,
	    n_e      = -1,
	    d        = -1,
	    n_edge   = -1,
	    n_n_edge = -1;

	const int p_p = GSL_MAX(1,p);
	plotting_element_info(&rst,&connect,&vtk_types,&connect_e,&n_n,&n_e,&d,&n_edge,&n_n_edge,p_p,e_type);

	struct Plotting_Nodes* p_nodes = calloc(1,sizeof *p_nodes); // returned
	struct Nodes* nodes = (struct Nodes*) p_nodes;

	nodes->rst = constructor_move_Matrix_d_d('C',n_n,d,true,rst); // destructed

	struct Vector_i** data_conn = malloc(n_e * sizeof *data_conn); // keep
	for (int i = 0; i < n_e; ++i) {
		const int n_conn = get_vtk_n_corners(vtk_types[i]);
		data_conn[i] = constructor_empty_Vector_i(n_conn); // keep
		for (int j = 0; j < n_conn; ++j)
			data_conn[i]->data[j] = connect[i*8+j];
	}
	free(connect);

	ptrdiff_t* exts_conn = malloc(1 * sizeof *exts_conn); // keep
	exts_conn[0] = n_e;
	p_nodes->connect = constructor_move_Multiarray_Vector_i_dyn_extents(1,exts_conn,true,data_conn); // destructed

	struct Vector_i** data_conn_e = malloc(n_edge * sizeof *data_conn_e); // keep
	for (int i = 0; i < n_edge; ++i) {
		data_conn_e[i] = constructor_empty_Vector_i(n_n_edge); // keep
		for (int j = 0; j < n_n_edge; ++j)
			data_conn_e[i]->data[j] = connect_e[i*n_n_edge+j];
	}
	free(connect_e);

	ptrdiff_t* exts_conn_e = malloc(1 * sizeof *exts_conn_e); // keep
	exts_conn_e[0] = n_edge;
	p_nodes->connect_e = constructor_move_Multiarray_Vector_i_dyn_extents(1,exts_conn_e,true,data_conn_e); // destructed

	p_nodes->vtk_types = constructor_move_Vector_i_i(n_e,true,vtk_types); // destructed

	p_nodes->vtk_types_e = constructor_empty_Vector_i(n_edge); // destructed
	set_to_value_Vector_i(p_nodes->vtk_types_e,VTK_POLY_LINE);

	return (const struct const_Plotting_Nodes*) p_nodes;
}

void destructor_const_Plotting_Nodes (const struct const_Plotting_Nodes*const p_nodes)
{
	destructor_Plotting_Nodes((struct Plotting_Nodes*)p_nodes);
}

void destructor_const_Plotting_Nodes_part (const struct const_Plotting_Nodes*const p_nodes)
{
	destructor_Plotting_Nodes_part((struct Plotting_Nodes*)p_nodes);
}

// Helper functions ************************************************************************************************* //

void print_Plotting_Nodes_tol (const struct Plotting_Nodes*const p_nodes, const double tol)
{
	struct Nodes* nodes = (struct Nodes*) p_nodes;

	printf("rst:\n");
	print_Matrix_d_tol(nodes->rst,tol);

	printf("connect:\n");
	print_Multiarray_Vector_i(p_nodes->connect);

	printf("vtk_types:\n");
	print_Vector_i(p_nodes->vtk_types);

	printf("vtk_types_e:\n");
	print_Vector_i(p_nodes->vtk_types_e);
}

void print_Plotting_Nodes (const struct Plotting_Nodes*const p_nodes)
{
	print_Plotting_Nodes_tol(p_nodes,EPS);
}

void print_const_Plotting_Nodes_tol (const struct const_Plotting_Nodes*const p_nodes, const double tol)
{
	print_Plotting_Nodes_tol((const struct Plotting_Nodes*const)p_nodes,tol);
}

void print_const_Plotting_Nodes (const struct const_Plotting_Nodes*const p_nodes)
{
	print_Plotting_Nodes((const struct Plotting_Nodes*const)p_nodes);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Get the number of edges of the input element type.
 *  \return See brief. */
static int get_n_e_type
	(const int e_type ///< \ref Element::type.
	);

static void destructor_Plotting_Nodes (struct Plotting_Nodes* p_nodes)
{
	destructor_Plotting_Nodes_part(p_nodes);

	struct Nodes* nodes = (struct Nodes*) p_nodes;

	assert(nodes->has_weights == false);
	destructor_Nodes(nodes);
}

static void destructor_Plotting_Nodes_part (struct Plotting_Nodes* p_nodes)
{
	destructor_Multiarray_Vector_i(p_nodes->connect);
	destructor_Multiarray_Vector_i(p_nodes->connect_e);
	destructor_Vector_i(p_nodes->vtk_types);
	destructor_Vector_i(p_nodes->vtk_types_e);
}

static void plotting_element_info
	(double **rst, int **connect, int **types, int **connect_e, int *n_n, int *n_e, int *d_, int *n_edge,
	 int *n_n_edge, const int p, const int e_type)
{
	const int P = p;
	int i, j, k, l, m, iMax, jMax, kMax, lMax, row, j0, jS, jC, count, N,
	             d, Nc, n_nOut, n_eOut, Ne,
	             layer, lBs, lTs, sum, u1 = 1,
	             *connectOut, *typesOut, *connect_eOut;
	double di, dj, dk, *rstOut;

	if (P == 0)
		printf("Error: Input P must be greater than 0 for plotting nodes.\n"), EXIT_MSG;

	Ne = get_n_e_type(e_type);
	const int NnEdge = ( e_type == LINE ? 1 : P+1 );

	connect_eOut = calloc(Ne*NnEdge , sizeof *connect_eOut); // keep

	if (e_type == LINE || e_type == QUAD || e_type == HEX) {
		int Indc,
		             N2,
		             nLINE[2], nQUAD[4], nHEX[8];
		int          sd, sP;
		double *r;

		// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
		d = 0;

		if      (e_type == LINE) d = 1;
		else if (e_type == QUAD) d = 2;
		else if (e_type == HEX)  d = 3;

		N = P+1;
		n_nOut = pow(N,d);
		n_eOut = pow(P,d);

		u1 = 1;
		sd = d;
		sP = P;

		r      = malloc(N       * sizeof *r);      // free
		rstOut = malloc(n_nOut*d * sizeof *rstOut); // keep

		connectOut  = calloc(n_eOut*8 , sizeof *connectOut); // keep
		typesOut    = malloc(n_eOut   * sizeof *typesOut);   // keep

		const struct const_Nodes* nodes = constructor_const_Nodes_tp(d,P,NODES_PLOT); // destructed
		for (i = 0; i < n_nOut*d; i++)
			rstOut[i] = nodes->rst->data[i];
		destructor_const_Nodes(nodes);

		free(r);

		N2 = pow(N,2);

		Indc = 0;
		for (k = 0, kMax = GSL_MAX(sP*GSL_MIN(sd-2,1),1); k < kMax; k++) {
			for (j = 0, jMax = GSL_MAX(P*GSL_MIN(d-1,u1),u1); j < jMax; j++) {
				for (i = 0; i < P; i++) {
					nLINE[0] = i;
					nLINE[1] = i+1;

					for (l = 0; l < 2; l++)             nQUAD[l] = nLINE[l]+N*j;
					for (l = 2, m = 1; l < 4; l++, m--) nQUAD[l] = nLINE[m] + N*(j+1);

					for (l = 0; l < 4; l++)             nHEX[l] = nQUAD[l] + N2*k;
					for (l = 4, m = 0; l < 8; l++, m++) nHEX[l] = nQUAD[m] + N2*(k+1);

					for (l = 0, lMax = pow(2,d); l < lMax; l++)
						connectOut[Indc*8+l] = nHEX[l];
					Indc++;
				}
			}
		}

		if (e_type == LINE) {
			for (i = 0; i < n_eOut; i++)
				typesOut[i] = 3;
			for (i = 0; i < Ne; i++) {
				switch (i) {
					case 0: j0 = 0; break;
					case 1: j0 = P; break;
				}
				connect_eOut[i] = j0;
			}
		} else if (e_type == QUAD) {
			for (i = 0; i < n_eOut; i++)
				typesOut[i] = 9;
			for (i = 0; i < Ne; i++) {
				switch (i) {
					case 0: j0 = 0;   jS = N; break;
					case 1: j0 = P;   jS = N; break;
					case 2: j0 = 0;   jS = 1; break;
					case 3: j0 = P*N; jS = 1; break;
					default:
						printf("Error: Unsupported.\n"), EXIT_MSG; break;
				}
				for (j = j0, count = 0; count < N; j += jS, count++)
					connect_eOut[i*N+count] = j;
			}
		} else if (e_type == HEX) {
			for (i = 0; i < n_eOut; i++)
				typesOut[i] = 12;
			for (i = 0; i < Ne; i++) {
				switch (i) {
					case 0:  j0 = 0;         jS = 1;   break;
					case 1:  j0 = P*N;       jS = 1;   break;
					case 2:  j0 = P*N*N;     jS = 1;   break;
					case 3:  j0 = (N*N-1)*N; jS = 1;   break;
					case 4:  j0 = 0;         jS = N;   break;
					case 5:  j0 = P;         jS = N;   break;
					case 6:  j0 = P*N*N;     jS = N;   break;
					case 7:  j0 = P*N*N+P;   jS = N;   break;
					case 8:  j0 = 0;         jS = N*N; break;
					case 9:  j0 = P;         jS = N*N; break;
					case 10: j0 = P*N;       jS = N*N; break;
					case 11: j0 = P*N+P;     jS = N*N; break;
					default:
						printf("Error: Unsupported.\n"), EXIT_MSG; break;
				}
				for (j = j0, count = 0; count < N; j += jS, count++)
					connect_eOut[i*N+count] = j;
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connect_e = connect_eOut;
		*n_n       = n_nOut;
		*n_e       = n_eOut;
		*d_        = d;
		*n_edge    = Ne;
		*n_n_edge  = NnEdge;
	} else if (e_type == TRI) {
		d = 2;
		Nc = 3;
		n_nOut = 1.0/2.0*(P+1)*(P+2);
		n_eOut = 0;
		for (i = 0; i < P; i++)
			n_eOut += 2*i+1;

		int iStart;
		double rst_V[Nc*d], BCoords[n_nOut*Nc];

		connectOut = calloc(n_eOut*8 , sizeof *connectOut); // keep
		typesOut   = malloc(n_eOut   * sizeof *typesOut); // keep

		// Determine barycentric coordinates of equally spaced nodes
		row = 0;
		for (j = 0; j <= P; j++) {
		for (i = 0, iMax = P-j; i <= iMax; i++) {
			di = i;
			dj = j;

			BCoords[row*Nc+2] = dj/P;
			BCoords[row*Nc+1] = di/P;
			BCoords[row*Nc+0] = 1.0 - (BCoords[row*Nc+1]+BCoords[row*Nc+2]);

			row++;
		}}

		// TRI vertex nodes
		rst_V[0*d+0] = -1.0; rst_V[0*d+1] = -1.0/sqrt(3.0);
		rst_V[1*d+0] =  1.0; rst_V[1*d+1] = -1.0/sqrt(3.0);
		rst_V[2*d+0] =  0.0; rst_V[2*d+1] =  2.0/sqrt(3.0);

		const struct const_Matrix_d* BCoords_M = constructor_move_const_Matrix_d_d('R',n_nOut,Nc,false,BCoords); // destructed
		const struct const_Matrix_d* rst_V_M   = constructor_move_const_Matrix_d_d('R',Nc,d,false,rst_V); // destructed
		struct Matrix_d* rst_M = constructor_mm_Matrix_d('N','N',1.0,BCoords_M,rst_V_M,'C'); // destructed
		rstOut = rst_M->data; // keep
		rst_M->owns_data = false;
		destructor_Matrix_d(rst_M);
		destructor_const_Matrix_d(BCoords_M);
		destructor_const_Matrix_d(rst_V_M);

		row = 0;
		// Regular TRIs
		for (j = P; j; j--) {
			iStart = 0;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+1;
				connectOut[row*8+2] = i+1+j;
				row++;
			}
		}

		// Inverted  TRIs
		for (j = P; j; j--) {
			iStart = 1;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j-1; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+j;
				connectOut[row*8+2] = i+j+1;
				row++;
			}
		}

		for (i = 0; i < n_eOut; i++)
			typesOut[i] = 5;

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = P; jS = P; break;
				case 1:  j0 = 0; jS = N; break;
				case 2:  j0 = 0; jS = 1; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connect_eOut[i*N+count] = j;

				if (j != j0 && jS > 1)
					jS--;
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connect_e = connect_eOut;
		*n_n       = n_nOut;
		*n_e       = n_eOut;
		*d_        = d;
		*n_edge    = Ne;
		*n_n_edge  = NnEdge;
	} else if (e_type == TET) {
		d = 3;
		Nc = 4;
		n_nOut = 1.0/6.0*(P+1)*(P+2)*(P+3);
		n_eOut = 0;
		for (i = 1; i <= P; i++) {
			// Regular TETs
			for (j = 1; j <= i; j++)
				n_eOut += j;

			// PYRs
			n_eOut += i*(i-1);

			// Inverted TETs
			if (i > 2)
				for (j = 1; j <= (i-2); j++)
					n_eOut += j;
		}

		double BCoords[n_nOut*Nc], rst_V[Nc*d];

		connectOut = calloc(n_eOut*8 , sizeof *connectOut); // keep
		typesOut   = malloc(n_eOut   * sizeof *typesOut); // keep

		// Determine barycentric coordinates of equally spaced nodes
		row = 0;
		for (k = 0; k <= P; k++) {
		for (j = 0, jMax = P-k; j <= jMax; j++) {
		for (i = 0, iMax = P-(j+k); i <= iMax; i++) {
			di = i;
			dj = j;
			dk = k;

			BCoords[row*Nc+3] = dk/P;
			BCoords[row*Nc+2] = dj/P;
			BCoords[row*Nc+1] = di/P;
			BCoords[row*Nc+0] = 1.0 - (BCoords[row*Nc+1]+BCoords[row*Nc+2]+BCoords[row*Nc+3]);

			row++;
		}}}

		// TET vertex nodes
		rst_V[0*d+0] = -1.0; rst_V[0*d+1] = -1.0/sqrt(3.0); rst_V[0*d+2] = -1.0/sqrt(6);
		rst_V[1*d+0] =  1.0; rst_V[1*d+1] = -1.0/sqrt(3.0); rst_V[1*d+2] = -1.0/sqrt(6);
		rst_V[2*d+0] =  0.0; rst_V[2*d+1] =  2.0/sqrt(3.0); rst_V[2*d+2] = -1.0/sqrt(6);
		rst_V[3*d+0] =  0.0; rst_V[3*d+1] =  0.0          ; rst_V[3*d+2] =  3.0/sqrt(6);

		const struct const_Matrix_d* BCoords_M = constructor_move_const_Matrix_d_d('R',n_nOut,Nc,false,BCoords); // destructed
		const struct const_Matrix_d* rst_V_M   = constructor_move_const_Matrix_d_d('R',Nc,d,false,rst_V); // destructed
		struct Matrix_d* rst_M = constructor_mm_Matrix_d('N','N',1.0,BCoords_M,rst_V_M,'C'); // destructed
		rstOut = rst_M->data; // keep
		rst_M->owns_data = false;
		destructor_Matrix_d(rst_M);
		destructor_const_Matrix_d(BCoords_M);
		destructor_const_Matrix_d(rst_V_M);

		row = 0;
		for (layer = P; layer; layer--) {
			lBs = 0;
			for (i = P; i > layer; i--) {
				sum = 0;
				for (j = 1; j <= i+1; j++)
					sum += j;
				lBs += sum;
			}

			lTs = 0;
			for (i = P; i > (layer-1); i--) {
				sum = 0;
				for (j = 1; j <= i+1; j++)
					sum += j;
				lTs += sum;
			}

			// Regular TETs
			for (j = 1; j <= layer; j++) {
			for (i = 0, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+0] = lBs + i + sum;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+0] + layer+2-j;

				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+3] = lTs + i + sum;

				row++;
			}}

			// PYRs
			for (j = 1, jMax = layer-1; j <= jMax; j++) {
			for (i = 1, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+4] = lBs + i + sum;
				connectOut[row*8+0] = connectOut[row*8+4] + layer+1-j;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;

				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+3] = lTs + i-1 + sum;
				connectOut[row*8+2] = connectOut[row*8+3] + 1;

				row++;

				connectOut[row*8+0] = connectOut[(row-1)*8+0];
				connectOut[row*8+1] = connectOut[(row-1)*8+1];
				connectOut[row*8+2] = connectOut[(row-1)*8+2];
				connectOut[row*8+3] = connectOut[(row-1)*8+3];
				connectOut[row*8+4] = connectOut[row*8+3] + layer+1-j;

				row++;
			}}

			// Inverted TETs
			for (j = 2, jMax = layer-1; j <= jMax; j++) {
			for (i = 1, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+0] = lBs + i +sum;

				sum = 0;
				for (k = 1, kMax = j-2; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+1] = lTs + i + sum;
				connectOut[row*8+2] = connectOut[row*8+1] + layer+1-j;
				connectOut[row*8+3] = connectOut[row*8+2] + 1;

				row++;
			}}
		}

		row = 0;
		for (i = P; i ; i-- ) {
			// Regular TETs
			sum = 0;
			for (j = 1; j <= i; j++)
				sum += j;
			for (jMax = sum; jMax--; ) {
				typesOut[row] = 10;
				row++;
			}

			// PYRs
			for (jMax = i*(i-1); jMax--; ) {
				typesOut[row] = 14;
				row++;
			}

			// Inverted TETs
			if (i > 2) {
				sum = 0;
				for (j = 1; j <= (i-2); j++)
					sum += j;

				for (jMax = sum; jMax--; ) {
					typesOut[row] = 10;
					row++;
				}
			}
		}

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = 0;             jS = 1;             jC = 1; break;
				case 1:  j0 = 0;             jS = N;             jC = 1; break;
				case 2:  j0 = 0;             jS = N*(N+1)/2.0;   jC = N; break;
				case 3:  j0 = P;             jS = P;             jC = 1; break;
				case 4:  j0 = P;             jS = N*(N+1)/2.0-1; jC = N; break;
				case 5:  j0 = N*(N+1)/2.0-1; jS = P*N/2.0;       jC = P; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connect_eOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC--;
				}
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connect_e = connect_eOut;
		*n_n       = n_nOut;
		*n_e       = n_eOut;
		*d_        = d;
		*n_edge    = Ne;
		*n_n_edge  = NnEdge;
	} else if (e_type == WEDGE) {
		d = 3;

		n_eOut = 0;
		for (i = 0; i < P; i++)
			n_eOut += 2*i+1;
		n_eOut *= P;

		double *rst_TRI, *rst_LINE;
		int n_n_TRI, n_n_LINE, NE_TRI, NE_LINE, tmp_i,
		             *connect_TRI, *connect_LINE, *dummy_types, *dummy_connect_e;

		plotting_element_info(&rst_TRI,&connect_TRI,&dummy_types,&dummy_connect_e,&n_n_TRI,&NE_TRI,&tmp_i,&tmp_i,&tmp_i,P,TRI); // free
		free(dummy_types); free(dummy_connect_e);
		plotting_element_info(&rst_LINE,&connect_LINE,&dummy_types,&dummy_connect_e,&n_n_LINE,&NE_LINE,&tmp_i,&tmp_i,&tmp_i,P,LINE); // free
		free(dummy_types); free(dummy_connect_e);

		n_nOut = n_n_TRI*n_n_LINE;

		rstOut     = malloc(n_nOut*d * sizeof *rstOut);     // keep
		connectOut = calloc(n_eOut*8 , sizeof *connectOut); // keep
		typesOut   = malloc(n_eOut   * sizeof *typesOut);   // keep

		row = 0;
		for (j = 0; j < n_n_LINE; j++) {
		for (i = 0; i < n_n_TRI; i++) {
			for (k = 0; k < 2; k++)
				rstOut[k*n_nOut+row] = rst_TRI[k*n_n_TRI+i];

			rstOut[2*n_nOut+row] = rst_LINE[j];

			row++;
		}}

		row = 0;
		for (j = 0, jMax = P; j < jMax; j++) {
		for (i = 0, iMax = n_eOut/P; i < iMax; i++) {
			for (k = 0; k < 3; k++)
				connectOut[row*8+k] = connect_TRI[i*8+k] + n_n_TRI*j;
			for (k = 0; k < 3; k++)
				connectOut[row*8+3+k] = connect_TRI[i*8+k] + n_n_TRI*(j+1);

			row++;
		}}

		for (i = 0; i < n_eOut; i++)
			typesOut[i] = 13;

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = P;               jS = P;  jC = 1; break;
				case 1:  j0 = 0;               jS = N;  jC = 1; break;
				case 2:  j0 = 0;               jS = 1;  jC = 1; break;
				case 3:  j0 = N*(N+1)/2.0*P+P; jS = P;  jC = 1; break;
				case 4:  j0 = N*(N+1)/2.0*P;   jS = N;  jC = 1; break;
				case 5:  j0 = N*(N+1)/2.0*P;   jS = 1;  jC = 1; break;
				case 6:  j0 = 0;               jS = N*(N+1)/2.0; jC = 0; break;
				case 7:  j0 = 1;               jS = N*(N+1)/2.0; jC = 0; break;
				case 8:  j0 = 2;               jS = N*(N+1)/2.0; jC = 0; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connect_eOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC--;
				}
			}
		}

		free(rst_TRI);
		free(connect_TRI);

		free(rst_LINE);
		free(connect_LINE);

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connect_e = connect_eOut;
		*n_n       = n_nOut;
		*n_e       = n_eOut;
		*d_        = d;
		*n_edge    = Ne;
		*n_n_edge  = NnEdge;
	} else if (e_type == PYR) {
		d = 3;
		Nc = 5;

		n_nOut = 0;
		for (i = 1, iMax = P+1; i <= iMax; i++)
			n_nOut += pow(i,2);

		n_eOut = 0;
		for (i = P; i ; i--) {
			n_eOut += pow(i,2);
			n_eOut += 2*i*(i-1);
			n_eOut += pow(i-1,2);
		}

		double *rst_QUAD;
		int n_n_QUAD, NE_QUAD, tmp_i,
		             *connect_QUAD, *dummy_types, *dummy_connect_e;

		rstOut     = malloc(n_nOut*d * sizeof *rstOut);     // keep
		connectOut = calloc(n_eOut*8 , sizeof *connectOut); // keep
		typesOut   = malloc(n_eOut   * sizeof *typesOut);   // keep

		row = 0;
		for (i = P; i ; i--) {
			di = i;
			plotting_element_info(&rst_QUAD,&connect_QUAD,&dummy_types,&dummy_connect_e,&n_n_QUAD,&NE_QUAD,&tmp_i,&tmp_i,&tmp_i,i,QUAD); // free
			free(dummy_types); free(dummy_connect_e);

			for (j = 0; j < n_n_QUAD; j++) {
				for (k = 0; k < 2; k++)
					rstOut[k*n_nOut+row] = rst_QUAD[k*n_n_QUAD+j]*di/P;

				rstOut[2*n_nOut+row] = -1.0/5.0*sqrt(2.0) + sqrt(2.0)/2.0*(1.0+(-1.0*di/P + 1.0*(P-di)/P));
				row++;
			}

			free(rst_QUAD);
			free(connect_QUAD);
		}
		for (k = 0; k < 2; k++)
			rstOut[k*n_nOut+row] = 0.0;
		rstOut[2*n_nOut+row] = 4.0/5.0*sqrt(2.0);

		row = 0;
		for (layer = P; layer; layer--) {
			lBs = 0;
			for (i = P; i > layer; i--)
				lBs += pow(i+1,2);

			lTs = 0;
			for (i = P; i > (layer-1); i--)
				lTs += pow(i+1,2);

			// Regular PYRs
			for (j = 0; j < layer; j++) {
			for (i = 0; i < layer; i++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+1] + layer+1;
				connectOut[row*8+3] = connectOut[row*8+2] - 1;
				connectOut[row*8+4] = lTs + i + j*layer;
				row++;
			}}

			// TETs
			for (j = 0; j < layer; j++) {
			for (i = 1; i < layer; i++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + layer+1;
				connectOut[row*8+2] = lTs + i-1 + j*layer;
				connectOut[row*8+3] = connectOut[row*8+2] + 1;
				row++;
			}}

			for (i = 0; i < layer; i++) {
			for (j = 1; j < layer; j++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = lTs + i + (j-1)*layer;
				connectOut[row*8+3] = connectOut[row*8+2] + layer;
				row++;
			}}

			// Inverted PYRs
			for (j = 1; j < layer; j++) {
			for (i = 1; i < layer; i++) {
				connectOut[row*8+0] = lTs + i-1 + (j-1)*layer;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+1] + layer;
				connectOut[row*8+3] = connectOut[row*8+2] - 1;
				connectOut[row*8+4] = lBs + i + j*(layer+1);
				row++;
			}}
		}

		row = 0;
		for (i = P; i ; i--) {
			// Regular PYRs
			for (j = 0, jMax = pow(i,2); j < jMax; j++) {
				typesOut[row] = 14;
				row++;
			}

			// TETs
			for (j = 0, jMax = 2*i*(i-1); j < jMax; j++) {
				typesOut[row] = 10;
				row++;
			}

			// Inverted PYRs
			for (j = 0, jMax = pow(i-1,2); j < jMax; j++) {
				typesOut[row] = 14;
				row++;
			}
		}

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = 0;     jS = 1;     jC = 0; break;
				case 1:  j0 = P*N;   jS = 1;     jC = 0; break;
				case 2:  j0 = 0;     jS = N;     jC = 0; break;
				case 3:  j0 = P;     jS = N;     jC = 0; break;
				case 4:  j0 = 0;     jS = N*N;   jC = 2*N-1; break;
				case 5:  j0 = P;     jS = N*N-1; jC = 2*N-1; break;
				case 6:  j0 = P*N;   jS = P*P+1; jC = 2*P-1; break;
				case 7:  j0 = N*N-1; jS = P*P;   jC = 2*P-1; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connect_eOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC -= 2;
				}
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connect_e = connect_eOut;
		*n_n       = n_nOut;
		*n_e       = n_eOut;
		*d_        = d;
		*n_edge    = Ne;
		*n_n_edge  = NnEdge;
	}
}

static int get_vtk_n_corners (const int vtk_type)
{
	switch (vtk_type) {
		case VTK_LINE:       return 2; break;
		case VTK_TRIANGLE:   return 3; break;
		case VTK_QUAD:       return 4; break;
		case VTK_TETRA:      return 4; break;
		case VTK_HEXAHEDRON: return 8; break;
		case VTK_WEDGE:      return 6; break;
		case VTK_PYR:        return 5; break;
		default:
			EXIT_ERROR("Error: Unsupported VTK_type (%d).\n",vtk_type);
			break;
	}
}

// Level 1 ********************************************************************************************************** //

static int get_n_e_type (const int e_type)
{
	switch (e_type) {
		case LINE:  return 2;  break;
		case TRI:   return 3;  break;
		case QUAD:  return 4;  break;
		case TET:   return 6;  break;
		case HEX:   return 12; break;
		case WEDGE: return 9;  break;
		case PYR:   return 8;  break;
		default: EXIT_ERROR("Unsupported: %d\n",e_type); break;
	}
}
