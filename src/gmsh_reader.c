#include <stdlib.h>
#include <stdio.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

#include "metis.h"
#include "parmetis.h"
#include "petscsys.h"
#include "mkl.h"

/*
 *	Purpose:
 *		Read in mesh file in gmsh format.
 *
 *	Comments:
 *
 *
 *		Periodic:
 *			The periodic processing is done in this function (in addition to simply reading the mesh file) because the
 *          desired output (PVe) cannot simply be read from the file (although this would be the ideal case).
 *
 *			Convention: *** IMPORTANT ***
 *
 *			For the code to properly recognize and match corresponding periodic vertices, it is necessary to follow the
 *			format specified below for the definition of Points, Lines, Surfaces, and Volumes in gmsh. This is not
 *			required when periodicity is not present.
 *
 *			Numbering (Geometry) :
 *
 *				Points   : All             1    - 1000
 *				Lines    : x (constant yz) 1001 - 2000
 *				           y (constant xz) 2001 - 3000
 *				           z (constant xy) 3001 - 4000
 *				Surfaces : xy (constant z) 4001 - 5000
 *				           xz (constant y) 5001 - 6000
 *				           yz (constant x) 6001 - 7000
 *				Volumes  : xyz             7001 - 8000
 *
 *			Numbering (Physical) :
 *				*0051 : Periodic (-x)
 *				*0052 : Periodic (+x)
 *				*0053 : Periodic (-y)
 *				*0054 : Periodic (+y)
 *				*0055 : Periodic (-z)
 *				*0056 : Periodic (+z)
 *
 *				* : 1 (Straight), 2 (Curved). This is the same as for all BCs
 *
 *			Periodic line and surface connections must also be specified.
 *
 *
 *		Parmetis:
 *
 *			Parmetis allows for weight constraints for mesh partitioning. When running using the full NS equations,
 *			check if more heavily weighting viscous elements results in improved performance; this would likely not be
 *			done in the initial mesh distribution but in a redistribution while the code is running. (ToBeDeleted)
 *
 *			ParMETIS_V3_PartMeshKway: While the manual indicates that part stores the local vertices, in this case, it
 *			                          stores the local elements (VOLUMEs).
 *
 *			Parmetis - 6) Restrictions & Limitations:
 *				Restriction 2 is false for the functions which are called here (i.e. Parmetis calls can be made using a
 *				single processor). There is, in fact, a special condition in the function calling metis if a single
 *				processor is used.
 *
 *
 *		Miscellaneous:
 *
 *			ToBeDeleted: After profiling the code, investigate using ParMETIS_V3_PartGeomKway if the current routine
 *			             turns out to be slow.
 *
 *			ToBeDeleted: The professor mentioned that the way OpenMPI is implemented here is "blocking", i.e. requires a
 *			             potentially large amount of wait time because of MPI_Barrier. There is a non-blocking way of
 *			             using OpenMPI sending the receiver information first. => Investigate later.
 *
 *	Notation:
 *		NVe     : (N)umber of (Ve)rtices
 *		NPVe    : (N)umber of (P)eriodic (Ve)rtices
 *		PVe     : List of (P)eriodic (Ve)rtices
 *		VeXYZ   : (Ve)rtex (xyz) locations
 *		NE      : (N)umber of (E)lements of each dimension (0D - 3D)
 *		NETotal : Sum of NE components.
 *
 *		EType  : Holds the (E)lement type numbering based on the gmsh convention as specified in Parameters.h
 *		ETags  : Holds the tags associated with each (E)lement.
 *		         1st tag: Physical tag (associated with all physical elements in gmsh) => Gives Boundary Condition info
 *		         2nd tag: Index of the geometry element in gmsh associated with this element (This is needed for the
 *		                  determination of periodic connections)
 *		EToVe  : (E)lement to (Ve)rtex correspondence.
 *		EToPrt : (E)lement to (P)a(rt)ition correspondence.
 *		         Note: The partitioning is done in parallel and the GLOBAL EToPrt array for all VOLUMEs is then
 *		               assembled on each processor using MPI.
 *
 *	References:
 *		Parmetis manual under: ${METIS_DIR}/manual
 *
 */

void gmsh_reader()
{
	// Initialize DB Parameters
	char  *MeshFile = DB.MeshFile;
	int    d        = DB.d,
	       Testing  = DB.Testing;
	int    PrintTesting = 0;

	// Standard datatypes
	char   StringRead[STRLEN_MAX], *strings, *stringe;
	int    i, j, k, iMax, dim, count, prt, flag, IndV, IndE, IndEV, IndP, ntags, type, *Nve, *Ed,
	       SectionNodes, SectionElements,
	       NVe, NETotal, *NE, *EType, *ETags, *EToVe, *EToPrt,
	       EPerProc, MPIsize, NVeRed, MPIrank,
	       NPeEnt, NPVe, PVeMax, *PVeOver, *PVeOverTrue, *PVe, *veMatch, PeEnt, NEnt, TagS, TagM, Ent, nodes[2],
	       *IndicesDummy, *NodesSOver, *NodesMOver, *NodesS, *NodesM, IndS, IndM, IndES, IndEM, E, NveMax,
	       NnS, NnM, *IndicesS, *IndicesM, *NodesSswap, *NodesMswap, IndPVe, Match,
	       Vs, *PVePossibleOver, *PVePossible, IndUnique, dimEntered[2], Es[4], IndVeXYZ[2];
	long   tmpl;
	double *VeXYZ, tmpd, *VeS, *VeM;

	FILE   *fID;

	struct S_ELEMENT *ELEMENT;

	// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
	NVe = -1; NPVe = 0; NETotal = -1; NveMax = -1;

	// Parmetis datatypes
	idx_t  *elmdist, *eptr, *eind, *numflag, *ncommonnodes,
	       *elmwgt = NULL,
		   *wgtflag, *ncon, *nparts, *options, *edgecut, *part;
	real_t *tpwgts, *ubvec;

	// MPI datatypes
	MPI_Status status;
	MPI_Comm comm = MPI_COMM_WORLD;


	NE = malloc(4 * sizeof *NE); //keep

	if ((fID = fopen(MeshFile,"r")) == NULL)
		printf("Mesh file: %s not present.\n",MeshFile), exit(1);

	// Find NVe, NETotal
	fscanf(fID,"%[^\n]\n",StringRead);
	while (!feof(fID)) {
		if (strstr(StringRead,"$Nodes") != NULL) {
			fscanf(fID,"%[^\n]\n",StringRead);
			sscanf(StringRead,"%d",&NVe);
		}

		if (strstr(StringRead,"$Elements") != NULL) {
			fscanf(fID,"%[^\n]\n",StringRead);
			sscanf(StringRead,"%d",&NETotal);
		}

		fscanf(fID,"%[^\n]\n",StringRead);
	}
	rewind(fID);

	VeXYZ = malloc(NVe*d     * sizeof *VeXYZ); // keep
	EType = malloc(NETotal*1 * sizeof *EType); // keep
	ETags = malloc(NETotal*2 * sizeof *ETags); // keep
	EToVe = malloc(NETotal*8 * sizeof *EToVe); // keep

	Nve = malloc(NETotal*1 * sizeof *Nve); // free
	Ed  = malloc(NETotal*1 * sizeof *Ed);  // free

	SectionNodes = 0;
	SectionElements = 0;

	fscanf(fID,"%[^\n]\n",StringRead);
	while (!feof(fID)) {
		fscanf(fID,"%[^\n]\n",StringRead);

		if (strstr(StringRead,"$Nodes") != NULL) {
			fscanf(fID,"%[^\n]\n",StringRead);
			fscanf(fID,"%[^\n]\n",StringRead);
			SectionNodes = 1;
		}
		if (strstr(StringRead,"$EndNodes") != NULL)
			SectionNodes = 0;

		if (SectionNodes) {
			strings = StringRead;

			// Vertex Index
			tmpl = strtol(strings,&stringe,10); strings = stringe;
			IndV = tmpl - 1;

			// Vertex Coordinates
			for (j = 0; j < d; j++) {
				tmpd = strtod(strings,&stringe); strings = stringe;
				VeXYZ[IndV*d+j] = tmpd;
			}
		}

		if (strstr(StringRead,"$Elements") != NULL) {
			fscanf(fID,"%[^\n]\n",StringRead);
			fscanf(fID,"%[^\n]\n",StringRead);
			SectionElements = 1;

			for (i = 0; i < 4; i++)
				NE[i] = 0;
		}
		if (strstr(StringRead,"$EndElements") != NULL)
			SectionElements = 0;

		// It is assumed that element types are grouped together in the msh file
		if (SectionElements) {
			strings = StringRead;

			// Element Index
			tmpl = strtol(strings,&stringe,10); strings = stringe;
			IndE = tmpl - 1;

			// Element Type
			EType[IndE] = strtol(strings,&stringe,10); strings = stringe;

			// Element Tags
			ntags = strtol(strings,&stringe,10); strings = stringe;
			if (ntags != 2)
				printf("Error: Too many node tags.\n"), exit(1);

			for (i = 0; i < ntags; i++) {
				tmpl = strtol(strings,&stringe,10); strings = stringe;
				ETags[IndE*ntags+i] = tmpl;
			}

			// Element Vertex Numbering
			ELEMENT = DB.ELEMENT; while(ELEMENT->type != EType[IndE]) ELEMENT = ELEMENT->next;
			type = ELEMENT->type;

			Nve[IndE] = ELEMENT->Nve;
			Ed[IndE]  = ELEMENT->d;   // Used in parmetis initialization

			for (i = 0; i < Nve[IndE]; i++) {
				tmpl = strtol(strings,&stringe,10); strings = stringe;
				EToVe[IndE*8+i] = tmpl - 1;
			}

			// Initialize remaining entries to NVe; needed for sorting (ToBeModified/Deleted)
			for (i = Nve[IndE]; i < 8; i++)
				EToVe[IndE*8+i] = NVe;

			// Number of each dimension
			if      (type == POINT)                                              NE[0]++;
			else if (type == LINE)                                               NE[1]++;
			else if (type == TRI || type == QUAD)                                NE[2]++;
			else if (type == TET || type == HEX || type == WEDGE || type == PYR) NE[3]++;
			else
				printf("Unsupported element type read from gmsh file.\n"), exit(1);
		}

		if (strstr(StringRead,"$Periodic") != NULL) {
			for (i = 0; i < 4; i++) Es[i] = 0;
			for (i = 0; i < 4; i++) {
			for (j = 0; j < i; j++) {
				Es[i] += NE[j];
			}}

			for (i = 0; i < 2; i++)
				dimEntered[i] = 0;

			fscanf(fID,"%[^\n]\n",StringRead);
			sscanf(StringRead,"%d",&NPeEnt);

			// Find Maximum number of possible periodic connections
			for (Vs = i = 0; i < d; i++)
				Vs += NE[i];

			PVePossibleOver = malloc(NETotal*pow(2,d-1) * sizeof *PVePossibleOver); // free
			for (i = 0, iMax = NETotal*pow(2,d-1); i < iMax; i++)
				PVePossibleOver[i] = NVe;

			for (E = 0, IndPVe = 0; E < Vs; E++) {
				if ((ETags[E*2+0] % 10000 >= 51) && (ETags[E*2+0] % 10000 <= 56)) {
					ELEMENT = DB.ELEMENT; while(ELEMENT->type != EType[E]) ELEMENT = ELEMENT->next;
					for (i = 0; i < ELEMENT->Nve; i++) {
						PVePossibleOver[IndPVe] = EToVe[E*8+i];
						IndPVe++;
					}
				}
			}
			PetscSortInt(IndPVe,&PVePossibleOver[0]);

			for (i = 1, IndUnique = 1; i < IndPVe; i++) {
				if (PVePossibleOver[i] != PVePossibleOver[IndUnique-1]) {
					PVePossibleOver[IndUnique] = PVePossibleOver[i];
					IndUnique++;
				}
			}
			free(PVePossibleOver);

			/*
			 * Slightly conservative estimate:
			 * Number of periodic vertices * Maximum number of periodic vertex connections +
			 * Maximum number of periodic corner vertices * Maximum number of corner periodic vertex connections
			 */
			PVeMax = (IndUnique*pow(2,d-1)+pow(2,d)*pow(2,d))/2;
			PVeOver = malloc(PVeMax*2 * sizeof *PVeOver); // free

			for (PeEnt = 0; PeEnt < NPeEnt; PeEnt++) {
				fscanf(fID,"%[^\n]\n",StringRead);

				sscanf(StringRead,"%d %d %d",&dim,&TagS,&TagM);
				if (dim == 0) { // Periodic 0D
					// gmsh provides the nodes => use directly

					fscanf(fID,"%[^\n]\n",StringRead);

NEnt = 10;
printf("%d\n\n\n",NEnt);


					sscanf(StringRead,"%d",&NEnt);

printf("%d\n\n\n",NEnt);

					for (Ent = 0; Ent < NEnt; Ent++) {
						fscanf(fID,"%[^\n]\n",StringRead);
						sscanf(StringRead,"%d %d",&nodes[0],&nodes[1]);
						for (i = 0; i < 2; i++) PVeOver[NPVe*2+i] = nodes[i]-1;
						PetscSortInt(2,&PVeOver[NPVe*2+0]);
						NPVe++;
					}
				} else if (dim == 1 || dim == 2) { // Periodic 1D or 2D
					// Find all possible periodic connections then eliminate connections unless d-1 coordinates are
					// shared.
					if ((dimEntered[0] == 0 && dim == 1 && d == 2) ||
					    (dimEntered[1] == 0 && dim == 2 && d == 3)) {

						find_periodic_connections(PVeOver,&NPVe,NVe);
//array_print_i(NPVe,2,PVeOver);

						PVeOverTrue = malloc(NPVe*2 * sizeof *PVeOverTrue); // free
						for (i = 0, iMax = NPVe*2; i < iMax; i++)
							PVeOverTrue[i] = PVeOver[i];

						veMatch = malloc(d * sizeof *veMatch); // free
						for (i = 0, IndP = 0; i < NPVe; i++) {
							for (j = 0, count = 0; j < d; j++) {
								veMatch[j] = 0;
								if (fabs(VeXYZ[PVeOver[i*2]*d+j] - VeXYZ[PVeOver[i*2+1]*d+j]) < NODETOL) {
									veMatch[j] = 1;
									count++;
								}
							}
							if (count == d-1) {
								for (j = 0; j < 2; j++)
									PVeOverTrue[IndP*2+j] = PVeOver[i*2+j];
								IndP++;
							}
						}
						free(veMatch);

//array_print_i(IndP,2,PVeOverTrue);
						for (i = IndP*2, iMax = NPVe*2; i < iMax; i++)
							PVeOver[i] = 0;

						NPVe = IndP;
						for (i = 0, iMax = NPVe*2; i < iMax; i++)
							PVeOver[i] = PVeOverTrue[i];
						free(PVeOverTrue);

//array_print_i(NPVe,2,PVeOver);
						dimEntered[dim-1] = 1;
					}

					IndicesDummy = malloc(2*NPVe * sizeof *IndicesDummy); // free
					array_sort_i(NPVe,2,PVeOver,IndicesDummy,'R','T');
					free(IndicesDummy);

					/*
					 * gmsh sometimes gives the nodes and sometimes does not... Further, even when gmsh lists nodes
					 * here, it does not necessarily list all of the nodes => find them and do not use the node list.
					 * (ToBeModified: Move to comments).
					 */
					if      (dim == 1) NveMax = 2;
					else if (dim == 2) NveMax = 4;

					NodesSOver = malloc(2*NE[dim]*NveMax * sizeof *NodesSOver); // free
					NodesMOver = malloc(2*NE[dim]*NveMax * sizeof *NodesMOver); // free
					for (i = 0; i < NE[dim]*NveMax; i++) {
						NodesSOver[i] = NVe;
						NodesMOver[i] = NVe;
					}

					// Make the list of slave and master nodes
					for (E = Es[dim], IndS = IndM = IndES = IndEM = 0; E < Es[dim+1]; E++) {
						if (ETags[E*2+1] == TagS) {
							for (j = 0; j < Nve[E]; j++) {
								NodesSOver[IndES*NveMax+j] = EToVe[E*8+j];
								IndS++;
							}
							IndES++;
						} else if (ETags[E*2+1] == TagM) {
							for (j = 0; j < Nve[E]; j++)  {
								NodesMOver[IndEM*NveMax+j] = EToVe[E*8+j];
								IndM++;
							}
							IndEM++;
						}
					}
					PetscSortInt(IndES*NveMax,&NodesSOver[0]);
					PetscSortInt(IndEM*NveMax,&NodesMOver[0]);

//array_print_i(1,IndS,NodesSOver);
//array_print_i(1,IndM,NodesMOver);

					// Remove duplicates
					for (i = 1, NnS = 1; i < IndS; i++) {
						if (NodesSOver[NnS-1] != NodesSOver[i]) {
							NodesSOver[NnS] = NodesSOver[i];
							NnS++;
						}
					}
					NodesS = malloc(NnS * sizeof *NodesS); // free
					for (i = 0; i < NnS; i++)
						NodesS[i] = NodesSOver[i];
					free(NodesSOver);

					for (i = 1, NnM = 1; i < IndM; i++) {
						if (NodesMOver[NnM-1] != NodesMOver[i]) {
							NodesMOver[NnM] = NodesMOver[i];
							NnM++;
						}
					}
					NodesM = malloc(NnM * sizeof *NodesM); // free
					for (i = 0; i < NnM; i++)
						NodesM[i] = NodesMOver[i];
					free(NodesMOver);

//array_print_i(1,NnS,NodesS);
//array_print_i(1,NnS,NodesM);

					if (NnS != NnM)
						printf("Error: NnS != NnM.\n"),	exit(1);

					// Coordinates which should match depending on which surface is periodic
					VeS = malloc(NnS*dim * sizeof *VeS); // free
					VeM = malloc(NnS*dim * sizeof *VeM); // free

					if (TagS > 1000 && TagS <= 2000) IndVeXYZ[0] = 0;                  // dim = 1: x-variable
					if (TagS > 2000 && TagS <= 3000) IndVeXYZ[0] = 1;                  // dim = 1: y-variable
					if (TagS > 3000 && TagS <= 4000) IndVeXYZ[0] = 2;                  // dim = 1: z-variable
					if (TagS > 4000 && TagS <= 5000) IndVeXYZ[0] = 0, IndVeXYZ[1] = 1; // dim = 2: xy-variable
					if (TagS > 5000 && TagS <= 6000) IndVeXYZ[0] = 0, IndVeXYZ[1] = 2; // dim = 2: xz-variable
					if (TagS > 6000 && TagS <= 7000) IndVeXYZ[0] = 1, IndVeXYZ[1] = 2; // dim = 2: yz-variable

					for (i = 0; i < NnS; i++) {
					for (j = 0; j < dim; j++) {
						VeS[i*dim+j] = VeXYZ[NodesS[i]*d+IndVeXYZ[j]];
						VeM[i*dim+j] = VeXYZ[NodesM[i]*d+IndVeXYZ[j]];
					}}

					IndicesS = malloc(NnS * sizeof *IndicesS); // free
					IndicesM = malloc(NnS * sizeof *IndicesM); // free
					for (i = 0; i < NnS; i++) {
						IndicesS[i] = i;
						IndicesM[i] = i;
					}
					array_sort_d(NnS,dim,VeS,IndicesS,'R','T');
					array_sort_d(NnS,dim,VeM,IndicesM,'R','T');

//array_print_d(NnS,dim,VeS);
//array_print_i(1,NnS,IndicesS);
//array_print_d(NnS,dim,VeM);
//array_print_i(1,NnS,IndicesM);

					free(VeS);
					free(VeM);

					NodesSswap = malloc(NnS * sizeof *NodesSswap); // free
					NodesMswap = malloc(NnS * sizeof *NodesMswap); // free
					for (i = 0; i < NnS; i++) {
						NodesSswap[i] = NodesS[i];
						NodesMswap[i] = NodesM[i];
					}

					for (i = 0; i < NnS; i++) {
						NodesS[i] = NodesSswap[IndicesS[i]];
						NodesM[i] = NodesMswap[IndicesM[i]];
					}
					free(NodesSswap);
					free(NodesMswap);
					free(IndicesS);
					free(IndicesM);

//array_print_i(1,NnS,NodesS);
//array_print_i(1,NnS,NodesM);

					PVePossible = malloc(NnS*2 * sizeof *PVePossible); // free
					for (i = 0; i < NnS; i++) {
						PVePossible[i*2+0] = NodesS[i];
						PVePossible[i*2+1] = NodesM[i];
						PetscSortInt(2,&PVePossible[i*2+0]);
					}
					free(NodesS);
					free(NodesM);

					IndicesDummy = malloc(NnS*2 * sizeof *IndicesDummy); // free
					array_sort_i(NnS,2,PVePossible,IndicesDummy,'R','T');
					free(IndicesDummy);

//array_print_i(NnS,2,PVePossible);

					for (i = 0, IndPVe = 0; i < NnS; i++) {
						while (IndPVe < NPVe && PVeOver[IndPVe*2+0] < PVePossible[i*2+0])
							IndPVe++;

						for (j = IndPVe, Match = 0; j < NPVe; j++) {
							if (PVeOver[j*2+0] == PVePossible[i*2+0] && PVeOver[j*2+1] == PVePossible[i*2+1])
								Match = 1;
						}

						if (!Match) {
							for (j = 0; j < 2; j++) PVeOver[NPVe*2+j] = PVePossible[i*2+j];
							NPVe++;
						}
					}
					free(PVePossible);

					// Jump to next periodic entity
					fscanf(fID,"%[^\n]\n",StringRead);
					sscanf(StringRead,"%d",&NEnt);
					for (Ent = 0; Ent < NEnt; Ent++) {
						fscanf(fID,"%[^\n]\n",StringRead);
					}
//array_print_i(NPVe,2,PVeOver);
				}
			}
			if (NPVe > PVeMax) {
				printf("Error: Too many periodic vertices found.\n");
				printf("NPVe: %d, PVeMax: %d\n",NPVe,PVeMax);
				exit(1);
			}

			// Store PVeOver in PVe and make links reflexive.
			PVe = malloc(2*NPVe*2 * sizeof *PVe); // keep
			for (k = 0; k < 2; k++) {
			for (i = 0; i < NPVe; i++) {
			for (j = 0; j < 2; j++) {
				if      (k == 0) PVe[k*NPVe*2+i*2+j] = PVeOver[i*2+j];
				else if (k == 1) PVe[k*NPVe*2+i*2+j] = PVeOver[i*2+1-j];
			}}}
			free(PVeOver);

			NPVe = 2*NPVe;
			IndicesDummy = malloc(NPVe*2 * sizeof *IndicesDummy); // free
			array_sort_i(NPVe,2,PVe,IndicesDummy,'R','T');
			free(IndicesDummy);
		}
	}
	if (NPVe == 0) PVe = malloc(0 * sizeof *PVe); // keep

	fclose(fID);

  printf("looking for linux bug...\n"); exit(1);

	// Add partition number to element list (in EToPrt) through communication between all processors
	MPIrank = DB.MPIrank;
	MPIsize = DB.MPIsize;

	// elmdist set based on gmsh element ordering
	EPerProc = (int) floor((double) NE[d]/MPIsize);

	elmdist = malloc((MPIsize+1) * sizeof *elmdist); // free

	elmdist[0] = 0;
	for (i = 1; i < MPIsize; i++ )
		elmdist[i] = elmdist[i-1] + EPerProc;
	elmdist[MPIsize] = NE[d];

	// count total number of vertices (including repeated) of all VOLUME elements
	eptr = malloc((NE[d]+1) * sizeof *eptr); // free

	// If elements are ALWAYS listed in order (by type), this can be simplified (ToBeModified)
	for (i = k = NVeRed = eptr[0] = 0; i < NETotal; i++) {
		if (Ed[i] == d) {
			NVeRed += Nve[i];
			eptr[k+1] = NVeRed;
			k++;
		}
	}

	eind = malloc(NVeRed * sizeof *eind); // free

	for (i = IndEV = 0; i < NETotal; i++) {
		if (Ed[i] == d) {
			for (k = eptr[IndEV], j = 0; k < eptr[IndEV+1]; k++) {
				eind[k] = EToVe[i*8+j];
				j++;
			}
			IndEV++;
		}
	}

	// Declare all variables on the heap so that the external function can see them
	numflag      = malloc(1 * sizeof *numflag); // free
	ncommonnodes = malloc(1 * sizeof *ncommonnodes); // free
	numflag[0]      = 0;
	ncommonnodes[0] = d;

//	elmwgt  = malloc(1 * sizeof *elmwgt); // free
	wgtflag = malloc(1 * sizeof *wgtflag); // free
	ncon    = malloc(1 * sizeof *ncon); // free
	nparts  = malloc(1 * sizeof *nparts); // free

//	elmwgt     = NULL;
	wgtflag[0] = 0;
	ncon[0]    = 1;
	nparts[0]  = MPIsize;

	tpwgts  = malloc(ncon[0]*MPIsize * sizeof *tpwgts); // free
	ubvec   = malloc(ncon[0]         * sizeof *ubvec); // free
	options = malloc(3               * sizeof *options); // free

	for (i = 0; i < ncon[0]; i++) {
		for (j = 0; j < MPIsize; j++) {
			tpwgts[i*MPIsize+j] = 1./((real_t) MPIsize);
		}
		ubvec[i] = 1.05; // Recommended in manual
	}

	for (i = 0; i < 3; i++)
		options[i] = 0; // See manual p.17 for options

	edgecut = malloc(1 *sizeof *edgecut); // free

	// Set part size according to number of elements initially placed on this proc
	part = malloc((elmdist[MPIrank+1]-elmdist[MPIrank]) * sizeof *part); // free

	ParMETIS_V3_PartMeshKway(
 		elmdist,eptr,eind,elmwgt,wgtflag,numflag,ncon,ncommonnodes,nparts,tpwgts,ubvec,options,edgecut,part,&comm);

	// Distribute partition information to all elements.

	// Initialize VOLUMEs on each proc
	EToPrt = malloc(NE[d] * sizeof *EToPrt); // keep
	flag = 1;

	for (i = 0; i < MPIsize; i++) {
		if (MPIrank == i) {
			count = 0;
			for (k = elmdist[MPIrank]; k < elmdist[MPIrank+1]; k++) {
				prt = part[count];

				for (j = 0; j < MPIsize; j++) {
					if (i != j) {
						MPI_Send(&k,1,MPI_INT,j,flag,comm);
						MPI_Send(&prt,1,MPI_INT,j,flag,comm);
					}
				}
				EToPrt[k] = prt;
				count++;
			}
			// bump partitions? This is taken from imex_adapt.c in Brian's code (ToBeModified)
			k = -1;
			prt = -1;
			for (j = 0; j < MPIsize; j++) {
				if (i != j) {
					MPI_Send(&k,1,MPI_INT,j,flag,comm);
					MPI_Send(&prt,1,MPI_INT,j,flag,comm);
				}
			}
		} else {
			k = 0;
			while (k >= 0) {
				MPI_Recv(&k,1,MPI_INT,i,flag,comm,&status);
				MPI_Recv(&prt,1,MPI_INT,i,flag,comm,&status);

				if (k >= 0) EToPrt[k] = prt;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Initialize imex type as explicit for all elements (ToBeModified)

	// Assign DB Parameters
	DB.NVe     = NVe;
	DB.NETotal = NETotal;
	DB.NE      = NE;
	DB.VeXYZ   = VeXYZ;
	DB.EType   = EType;
	DB.ETags   = ETags;
	DB.EToVe   = EToVe;
	DB.EToPrt  = EToPrt;
	DB.NPVe    = NPVe;
	DB.PVe     = PVe;

	// Testing
	if (!(MPIrank) && Testing && PrintTesting) {
		printf("NE:\n");     array_print_i(1,4,DB.NE);
		printf("VeXYZ:\n");  array_print_d(DB.NVe,d,DB.VeXYZ);
		printf("EType:\n");  array_print_i(1,DB.NETotal,DB.EType);
		printf("ETags:\n");  array_print_i(DB.NETotal,2,DB.ETags);
		printf("EToVe:\n");  array_print_i(DB.NETotal,8,DB.EToVe);
		printf("EToPrt:\n"); array_print_i(1,DB.NE[d],DB.EToPrt);
		printf("PVe:\n");    array_print_i(DB.NPVe,2,DB.PVe);
	}

	// Free memory
	free(Nve);
	free(Ed);

	METIS_Free(elmdist);
	METIS_Free(eptr);
	METIS_Free(eind);
//	METIS_Free(elmwgt);
	METIS_Free(wgtflag);
	METIS_Free(numflag);
	METIS_Free(ncon);
	METIS_Free(ncommonnodes);
	METIS_Free(nparts);
	METIS_Free(tpwgts);
	METIS_Free(ubvec);
	METIS_Free(options);
	METIS_Free(edgecut);
	METIS_Free(part);
}
