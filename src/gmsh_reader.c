// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "gmsh_reader.h"

#include <stdlib.h>
#include <stdio.h>

#include "metis.h"
#include "parmetis.h"
#include "petscsys.h"
#include "mkl.h"
 
#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"

#include "array_sort.h"
#include "array_print.h"
#include "array_find_index.h"

/*
 *	Purpose:
 *		Read in mesh file in gmsh format.
 *
 *	Comments:
 *
 *
 *		Periodic:
 *			The periodic processing is done in this function (in addition to simply reading the mesh file) because the
 *			desired output (PVe) cannot simply be read from the file (although this would be the ideal case).
 *
 *			Convention: *** IMPORTANT ***
 *
 *			For the code to properly recognize and match corresponding periodic vertices, it is necessary to follow the
 *			format specified below for the definition of Points, Lines, Surfaces, and Volumes in gmsh. This is not
 *			required when periodicity is not present. In this case, Lines and Surfaces are generally defined with
 *			arbitrary numbering within the appropriate range.
 *			Note: The numbering convention can be changed by changing the appropriate magic numbers in parameters.h.
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

void find_periodic_connections(unsigned int *PVe, unsigned int *NPVePointer, const unsigned int VeMax);

void gmsh_reader(void)
{
	// Initialize DB Parameters
	char         *MeshFile = DB.MeshFile;
	unsigned int d         = DB.d;

	unsigned int PrintTesting = 0;

	// Standard datatypes
	char         StringRead[STRLEN_MAX], *strings, *stringe;
	unsigned int i, j, k, iMax, dim, count, flag, IndV, IndE, IndEV, IndP, ntags, type, *Nve, *Ed, *VeCGmsh,
	             SectionNodes, SectionElements,
	             NVe, NETotal, *NE, *EType, *ETags, *EToVe, *EToPrt,
	             EPerProc, MPIsize, NVeRed, MPIrank,
	             NPeEnt, NPVe, PVeMax, *PVeOver, *PVeOverTrue, *PVe, *veMatch, PeEnt, NEnt, TagS, TagM, Ent, nodes[2],
	             *IndicesDummy, *NodesSOver, *NodesMOver, *NodesS, *NodesM, IndS, IndM, IndES, IndEM, E, NveMax,
	             NnS, NnM, *IndicesS, *IndicesM, *NodesSswap, *NodesMswap, IndPVe, Match,
	             Vs, *PVePossibleOver, *PVePossible, IndUnique, dimEntered[2], Es[4], IndVeXYZ[2];
	int          is, ks, prts;
	long         tmpl;
	double       *VeXYZ, tmpd, *VeS, *VeM;

	FILE         *fID;

	struct S_ELEMENT *ELEMENT;

	// silence
	NVe = 0; NPVe = 0; NETotal = 0; NveMax = 0;
	PVe = malloc(0 * sizeof *PVe); // keep (freed and reallocated if NPVe > 0 after reading file)
	IndVeXYZ[0] = 0;
	IndVeXYZ[1] = 0;

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
		printf("Mesh file: %s not present.\n",MeshFile), EXIT_MSG;

	// Find NVe, NETotal
	while (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		if (strstr(StringRead,"$Nodes")) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
				sscanf(StringRead,"%d",&NVe);
		}

		if (strstr(StringRead,"$Elements")) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
				sscanf(StringRead,"%d",&NETotal);
		}
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

	while (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		if (strstr(StringRead,"$Nodes")) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
			SectionNodes = 1;
		}
		if (strstr(StringRead,"$EndNodes"))
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

		if (strstr(StringRead,"$Elements")) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
				SectionElements = 1;

			for (i = 0; i < 4; i++)
				NE[i] = 0;
		}
		if (strstr(StringRead,"$EndElements"))
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
			VeCGmsh = ELEMENT->VeCGmsh;

			Nve[IndE] = ELEMENT->Nve;
			Ed[IndE]  = ELEMENT->d;   // Used in parmetis initialization

			for (i = 0; i < Nve[IndE]; i++) {
				tmpl = strtol(strings,&stringe,10); strings = stringe;
				EToVe[IndE*8+VeCGmsh[i]] = tmpl - 1;
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

		if (strstr(StringRead,"$Periodic")) {
			for (i = 0; i < 4; i++) Es[i] = 0;
			for (i = 0; i < 4; i++) {
			for (j = 0; j < i; j++) {
				Es[i] += NE[j];
			}}

			for (i = 0; i < 2; i++)
				dimEntered[i] = 0;

			if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
				sscanf(StringRead,"%d",&NPeEnt);

			// Find Maximum number of possible periodic connections
			for (Vs = i = 0; i < d; i++)
				Vs += NE[i];

			PVePossibleOver = malloc(NETotal*pow(2,d-1) * sizeof *PVePossibleOver); // free
			for (i = 0, iMax = NETotal*pow(2,d-1); i < iMax; i++)
				PVePossibleOver[i] = NVe;

			for (E = 0, IndPVe = 0; E < Vs; E++) {
				if ((ETags[E*2+0] % BC_STEP_SC >= PERIODIC_XL) && (ETags[E*2+0] % BC_STEP_SC <= PERIODIC_ZR)) {
					ELEMENT = DB.ELEMENT; while(ELEMENT->type != EType[E]) ELEMENT = ELEMENT->next;
					for (i = 0; i < ELEMENT->Nve; i++) {
						PVePossibleOver[IndPVe] = EToVe[E*8+i];
						IndPVe++;
					}
				}
			}
			PetscSortInt(IndPVe,(int *)PVePossibleOver);

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
				if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
					sscanf(StringRead,"%d %d %d",&dim,&TagS,&TagM);
				if (dim == 0) { // Periodic 0D
					// gmsh provides the nodes => use directly

					if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
					// skip line beginnning with "Affine" if present
					if (strstr(StringRead,"Affine"))
						if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }

					sscanf(StringRead,"%d",&NEnt);

					for (Ent = 0; Ent < NEnt; Ent++) {
						if (fscanf(fID,"%[^\n]\n",StringRead) == 1)
							sscanf(StringRead,"%d %d",&nodes[0],&nodes[1]);
						for (i = 0; i < 2; i++)
							PVeOver[NPVe*2+i] = nodes[i]-1;
						PetscSortInt(2,(int *)&PVeOver[NPVe*2+0]);
						NPVe++;
					}
				} else if (dim == 1 || dim == 2) { // Periodic 1D or 2D
					// Find all possible periodic connections then eliminate connections unless d-1 coordinates are
					// shared.
					if ((dimEntered[0] == 0 && dim == 1 && d == 2) ||
					    (dimEntered[1] == 0 && dim == 2 && d == 3)) {

						find_periodic_connections(PVeOver,&NPVe,NVe);
//array_print_i(NPVe,2,PVeOver,'R');

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

//array_print_i(IndP,2,PVeOverTrue,'R');
						for (i = IndP*2, iMax = NPVe*2; i < iMax; i++)
							PVeOver[i] = 0;

						NPVe = IndP;
						for (i = 0, iMax = NPVe*2; i < iMax; i++)
							PVeOver[i] = PVeOverTrue[i];
						free(PVeOverTrue);

//array_print_i(NPVe,2,PVeOver,'R');
						dimEntered[dim-1] = 1;
					}

					IndicesDummy = malloc(2*NPVe * sizeof *IndicesDummy); // free
					array_sort_ui(NPVe,2,PVeOver,IndicesDummy,'R','T');
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
					PetscSortInt(IndES*NveMax,(int *)NodesSOver);
					PetscSortInt(IndEM*NveMax,(int *)NodesMOver);

//array_print_i(1,IndS,NodesSOver,'R');
//array_print_i(1,IndM,NodesMOver,'R');

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

//array_print_i(1,NnS,NodesS,'R');
//array_print_i(1,NnS,NodesM,'R');

					if (NnS != NnM)
						printf("Error: NnS != NnM.\n"), exit(1);

					// Coordinates which should match depending on which surface is periodic
					VeS = malloc(NnS*dim * sizeof *VeS); // free
					VeM = malloc(NnS*dim * sizeof *VeM); // free

					if (TagS >= GMSH_XLINE_MIN  && TagS < GMSH_YLINE_MIN)  IndVeXYZ[0] = 0;
					if (TagS >= GMSH_YLINE_MIN  && TagS < GMSH_ZLINE_MIN)  IndVeXYZ[0] = 1;
					if (TagS >= GMSH_ZLINE_MIN  && TagS < GMSH_XYFACE_MIN) IndVeXYZ[0] = 2;
					if (TagS >= GMSH_XYFACE_MIN && TagS < GMSH_XZFACE_MIN) IndVeXYZ[0] = 0, IndVeXYZ[1] = 1;
					if (TagS >= GMSH_XZFACE_MIN && TagS < GMSH_YZFACE_MIN) IndVeXYZ[0] = 0, IndVeXYZ[1] = 2;
					if (TagS >= GMSH_YZFACE_MIN && TagS < GMSH_XYZVOL_MIN) IndVeXYZ[0] = 1, IndVeXYZ[1] = 2;

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

//array_print_d(NnS,dim,VeS,'R');
//array_print_i(1,NnS,IndicesS,'R');
//array_print_d(NnS,dim,VeM,'R');
//array_print_i(1,NnS,IndicesM,'R');

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

//array_print_i(1,NnS,NodesS,'R');
//array_print_i(1,NnS,NodesM,'R');

					PVePossible = malloc(NnS*2 * sizeof *PVePossible); // free
					for (i = 0; i < NnS; i++) {
						PVePossible[i*2+0] = NodesS[i];
						PVePossible[i*2+1] = NodesM[i];
						PetscSortInt(2,(int *)&PVePossible[i*2+0]);
					}
					free(NodesS);
					free(NodesM);

					IndicesDummy = malloc(NnS*2 * sizeof *IndicesDummy); // free
					array_sort_ui(NnS,2,PVePossible,IndicesDummy,'R','T');
					free(IndicesDummy);

//array_print_i(NnS,2,PVePossible,'R');

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
					if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
					// skip line beginnning with "Affine" if present
					if (strstr(StringRead,"Affine"))
						if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
					sscanf(StringRead,"%d",&NEnt);
					for (Ent = 0; Ent < NEnt; Ent++) {
						if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
					}
//array_print_i(NPVe,2,PVeOver,'R');
				}
			}
			if (NPVe > PVeMax) {
				printf("Error: Too many periodic vertices found.\n");
				printf("NPVe: %d, PVeMax: %d\n",NPVe,PVeMax);
				exit(1);
			}

			// Store PVeOver in PVe and make links reflexive.
			free(PVe);
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
			array_sort_ui(NPVe,2,PVe,IndicesDummy,'R','T');
			free(IndicesDummy);
		}
	}
	fclose(fID);

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
	for (i = 0, k = 0, NVeRed = 0, eptr[0] = 0; i < NETotal; i++) {
		if (Ed[i] == d) {
			NVeRed += Nve[i];
			eptr[k+1] = NVeRed;
			k++;
		}
	}

	eind = malloc(NVeRed * sizeof *eind); // free

	for (i = 0, IndEV = 0; i < NETotal; i++) {
		if (Ed[i] == d) {
			for (ks = eptr[IndEV], j = 0; ks < eptr[IndEV+1]; ks++) {
				eind[ks] = EToVe[i*8+j];
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

	for (is = 0; is < ncon[0]; is++) {
		for (j = 0; j < MPIsize; j++) {
			tpwgts[is*MPIsize+j] = 1./((real_t) MPIsize);
		}
		ubvec[is] = 1.05; // Recommended in manual
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
			for (ks = elmdist[MPIrank]; ks < elmdist[MPIrank+1]; ks++) {
				prts = part[count];

				for (j = 0; j < MPIsize; j++) {
					if (i != j) {
					//	printf("%d\n",MPI_Send(&ks,1,MPI_INT,j,flag,comm));
						if (MPI_Send(&ks,1,MPI_INT,j,flag,comm) != 1)
							printf("Error: MPI_Send error.\n"), exit(1);
						if (MPI_Send(&prts,1,MPI_INT,j,flag,comm) != 1)
							printf("Error: MPI_Send error.\n"), exit(1);
					}
				}
				EToPrt[ks] = prts;
				count++;
			}
			// bump partitions? This is taken from imex_adapt.c in Brian's code (ToBeModified)
			ks = -1;
			prts = -1;
			for (j = 0; j < MPIsize; j++) {
				if (i != j) {
					if (MPI_Send(&ks,1,MPI_INT,j,flag,comm) != 1)
						printf("Error: MPI_Send error.\n"), exit(1);
					if (MPI_Send(&prts,1,MPI_INT,j,flag,comm) != 1)
						printf("Error: MPI_Send error.\n"), exit(1);
				}
			}
		} else {
			ks = 0;
			while (ks >= 0) {
				if (MPI_Recv(&ks,1,MPI_INT,i,flag,comm,&status) != 1)
					printf("Error: MPI_Recv error.\n"), exit(1);
				if (MPI_Recv(&prts,1,MPI_INT,i,flag,comm,&status) != 1)
					printf("Error: MPI_Recv error.\n"), exit(1);

				if (ks >= 0) EToPrt[ks] = prts;
			}
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Initialize imex type as explicit for all elements (ToBeModified)

	// Assign DB Parameters
	DB.NVe     = NVe;
	DB.NPVe    = NPVe;
	DB.PVe     = PVe;
	DB.VeXYZ   = VeXYZ;
	DB.NE      = NE;
	DB.NETotal = NETotal;
	DB.EType   = EType;
	DB.ETags   = ETags;
	DB.EToVe   = EToVe;
	DB.EToPrt  = EToPrt;

	if (!DB.MPIrank && PrintTesting) {
		printf("NE:\n");     array_print_ui(1,4,DB.NE,'R');
		printf("VeXYZ:\n");  array_print_d(DB.NVe,d,DB.VeXYZ,'R');
		printf("EType:\n");  array_print_ui(1,DB.NETotal,DB.EType,'R');
		printf("ETags:\n");  array_print_ui(DB.NETotal,2,DB.ETags,'R');
		printf("EToVe:\n");  array_print_ui(DB.NETotal,8,DB.EToVe,'R');
		printf("EToPrt:\n"); array_print_ui(1,DB.NE[d],DB.EToPrt,'R');
		printf("PVe:\n");    array_print_ui(DB.NPVe,2,DB.PVe,'R');
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

void find_periodic_connections(unsigned int *PVe, unsigned int *NPVePointer, const unsigned int VeMax)
{
	/*
	 *	Purpose:
	 *		Given a list of periodic vertex correspondence, find all possible matches between the vertices.
	 *
	 *	Comments:
	 *		Each row of the list is sorted in ascending order. Hence, reversed entries are redundant and not included.
	 *
	 *	Notation:
	 *		PVe : List of periodic vertices (pve x 2 array)
	 *		pve : Index of last row in PVe
	 *
	 *	References:
	 */

	unsigned int i, j, k,
	             pve, *IndicesDummy, Modified, *PVe1D, NUnique, *PVeUnique, *PVeUniqueOver, *PVeMatches, *IndPVeMatches,
	             IndRow, LenF, match, row, col, n2, *PVePotential, IndPVe;

	pve = *NPVePointer;

// array_print_i(pve,2,PVe,'R');

	// Sort Rows
	for (i = 0; i < pve; i++)
		PetscSortInt(2,(int *)&PVe[i*2+0]);

	// Sort Columns
	IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
	array_sort_ui(pve,2,PVe,IndicesDummy,'R','T');
	free(IndicesDummy);

// array_print_i(pve,2,PVe,'R');

	Modified = 1;
	while (Modified) {
		Modified = 0;

		PVe1D = malloc(pve*2 * sizeof *PVe1D); // free
		for (i = k = 0; i < pve; i++) {
		for (j = 0; j < 2; j++) {
			PVe1D[k] = PVe[i*2+j];
			k++;
		}}

		PetscSortInt(k,(int *)PVe1D);
// array_print_i(1,k,PVe1D,'R');

		PVeUniqueOver = malloc(pve*2 * sizeof *PVeUniqueOver); // free
		PVeUniqueOver[0] = PVe1D[0];
		for (i = 1, NUnique = 1; i < pve*2; i++) {
			if (PVe1D[i] != PVeUniqueOver[NUnique-1]) {
				PVeUniqueOver[NUnique] = PVe1D[i];
				NUnique++;
			}
		}
		free(PVe1D);

		PVeUnique = malloc(NUnique * sizeof *PVeUnique); // free
		for (i = 0; i < NUnique; i++)
			PVeUnique[i] = PVeUniqueOver[i];
		free(PVeUniqueOver);

// array_print_i(1,NUnique,PVeUnique,'R');

		PVeMatches    = malloc(NUnique*8 * sizeof *PVeMatches); // free
		IndPVeMatches = malloc(NUnique   * sizeof *IndPVeMatches); // free
		for (i = 0; i < NUnique*8; i++) PVeMatches[i]    = VeMax;
		for (i = 0; i < NUnique; i++)   IndPVeMatches[i] = 0;

		for (k = 0; k < 2; k++) {
		for (i = 0; i < pve; i++) {
			array_find_indexo_ui(NUnique,PVeUnique,PVe[i*2+k],&IndRow,&LenF);

			for (j = 0, match = 0; j < IndPVeMatches[IndRow]; j++) {
				if (PVe[i*2+1-k] == PVeMatches[IndRow*8+j])
					match = 1;
			}

			if (!match) {
				PVeMatches[IndRow*8+IndPVeMatches[IndRow]] = PVe[i*2+1-k];

				IndPVeMatches[IndRow]++;
				if (IndPVeMatches[IndRow] > 7)
					printf("Error: Too many periodic connections\n"), exit(1);
			}
		}}

		// sort rows
		for (i = 0; i < NUnique; i++)
			PetscSortInt(IndPVeMatches[i]+1,(int *)&PVeMatches[i*8+0]);

// array_print_i(NUnique,8,PVeMatches,'R');
// array_print_i(1,NUnique,IndPVeMatches,'R');

		PVePotential = malloc(2 * sizeof *PVePotential); // free

		for (row = 0; row < NUnique; row++) {
		for (col = 0; col < IndPVeMatches[row]; col++) {
			array_find_indexo_ui(NUnique,PVeUnique,PVeMatches[row*8+col],&IndRow,&LenF);
			for (i = 0; i < IndPVeMatches[IndRow]; i++) {
				n2 = PVeMatches[IndRow*8+i];
				if (PVeUnique[row] != n2) {
					PVePotential[0] = PVeUnique[row];
					PVePotential[1] = n2;
					PetscSortInt(2,(int *)PVePotential);

					for (IndPVe = 0, match = 0; IndPVe < pve; IndPVe++) {
						if (PVe[IndPVe*2+0] == PVePotential[0] &&
							PVe[IndPVe*2+1] == PVePotential[1]) {
								match = 1;
								break;
						}
					}

					if (!match) {
						for (j = 0; j < 2; j++)
							PVe[pve*2+j] = PVePotential[j];

						pve++;
						Modified = 1;
					}
				}
			}
		}}
		free(PVeUnique);
		free(PVePotential);

		IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
		array_sort_ui(pve,2,PVe,IndicesDummy,'R','T');
		free(IndicesDummy);

// array_print_i(pve,2,PVe,'R');

		free(PVeMatches);
		free(IndPVeMatches);
	}

	*NPVePointer = pve;
}
