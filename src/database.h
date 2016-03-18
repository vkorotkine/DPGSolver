/*
  Purpose:
    Define global parameters and objects
 
  Comments:
    The notation is presented in the first routine in which parameters appear.
 
  Notation:
 
  References:
 
*/

#ifndef DPG__database_h__INCLUDED
#define DPG__database_h__INCLUDED

struct S_DB {
  // MPI and PETSC
  int MPIsize, MPIrank;

  // Initialization
  char *prefix, *MeshFile, *MeshPath, *MeshType, *Form, *NodeType, *BasisType, *TestCase;
  int  d, ML, Vectorized, EFE, Collocated, Adaptive, P, PMax, Restart, Testing;

  // Parameters
  int  NP, NDE, **SF_BE, AC, ExactGeom,
       PR, PP, PGs, *PGc, **PCs, **PCc, **PJs, **PJc,
       *PF, **PFrs, **PFrc, **PIfs, **PIfc, **PIvs, **PIvc;
  char *Parametrization,
       ***NodeTypeS,   ***NodeTypeF,   ***NodeTypeFrs, ***NodeTypeFrc,
       ***NodeTypeIfs, ***NodeTypeIfc, ***NodeTypeIvs, ***NodeTypeIvc;

  // Mesh
  struct S_ELEMENT *ELEMENT;

  int    NVe, NPVe, NfMax, NfveMax, *PVe, NETotal, *NE, *EType, *ETags, *EToVe, *EToPrt,
         NV, NGF, NVC, NGFC, *VToV, *VToF, *VToGF, *VToBC, *GFToVe, *VC, *GFC;
  double *VeXYZ;

};
extern struct S_DB DB;

typedef struct S_ELEMENT {
  // Mesh
  int present, type, d, Nve, Nf, *Nfve, *VeC, *VeE, *VeF;
  //int present, type, d, Nve, Nf, Nfve[2], VeC[8], VeE[12*2], VeF[6*4];

  struct S_ELEMENT *next;
} S_ELEMENT;

#endif // DPG__database_h_INCLUDED
