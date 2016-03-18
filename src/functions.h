/*
  Purpose:
    Set function prototypes
 
  Comments:
 
  Notation:
 
  References:
 
*/

#ifndef DPG__functions_h__INCLUDED
#define DPG__functions_h__INCLUDED

// Preprocessor
extern void Initialization              (int nargc, char **argv);
extern void SetupParameters             (void);
extern void SetupMesh                   (void);
extern void   GmshReader                (void);
extern void   SetupConnectivity         (void);
extern void   SetupPeriodic             (void);
extern void     FindPeriodicConnections (int *Pve, int *pvePointer, int VeMax);
extern void SetupOperators              (void);
extern void   Cubature_TP               (double **xir, double **W, int **Con, int *Nn, int *ToReturn, int P, int d, char *NodeType);


// Extern structs
extern struct S_ELEMENT *New_ELEMENT (void);


// Memory Management
extern void MemoryFree         (void);
extern void MemoryConstructors (void);
extern void MemoryDestructorE  (struct S_ELEMENT *ELEMENT);


// Array Processing
  // Sorting
  extern void ArraySorti       (int NRows, int NCols, int    *A, int *Indices, char ordering, char trans);
  extern void ArraySortd       (int NRows, int NCols, double *A, int *Indices, char ordering, char trans);
  extern void ArrayFindIndexOi (int LenA, int *A, int val, int *IdxF, int *LenF);

  // Norms
  extern double ArrayNormd (int LenA, double *A, char *NormType);

  // Printing
  extern void ArrayPrinti  (int m, int n, int *A);
  extern void ArrayPrintl  (int m, int n, long *A);
  extern void ArrayPrintll (int m, int n, long long *A);
  extern void ArrayPrintf  (int m, int n, float *A);
  extern void ArrayPrintd  (int m, int n, double *A);
  extern void ArrayPrintld (int m, int n, long double *A);

  // Memory Management
  extern void ArrayFree2c  (int iMax, char **A);
  extern void ArrayFree2i  (int iMax, int **A);
  extern void ArrayFree2l  (int iMax, long **A);
  extern void ArrayFree2ll (int iMax, long long **A);
  extern void ArrayFree2f  (int iMax, float **A);
  extern void ArrayFree2d  (int iMax, double **A);
  extern void ArrayFree2ld (int iMax, long double **A);
  extern void ArrayFree3c  (int iMax, int jMax, char ***A);
  extern void ArrayFree3i  (int iMax, int jMax, int ***A);
  extern void ArrayFree3l  (int iMax, int jMax, long ***A);
  extern void ArrayFree3ll (int iMax, int jMax, long long ***A);
  extern void ArrayFree3f  (int iMax, int jMax, float ***A);
  extern void ArrayFree3d  (int iMax, int jMax, double ***A);
  extern void ArrayFree3ld (int iMax, int jMax, long double ***A);



#endif // DPG_functions_h__INCLUDED
