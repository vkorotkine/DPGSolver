// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__bases_h__INCLUDED
#define DPG__bases_h__INCLUDED

#include "matrix_structs.h"

typedef double *(*basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                               const unsigned int d);
typedef double **(*grad_basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn,
                                     unsigned int *NbfOut, const unsigned int d);

extern double jacobiP          (const double x, const double alpha, const double beta, const int N);
extern double grad_jacobiP     (const double x, const double alpha, const double beta, const int N);

extern double *basis_TP        (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double *basis_SI        (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double *basis_PYR       (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_TP  (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_SI  (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_PYR (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern void   rst_to_abc_SI    (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);
extern void   rst_to_abc_PYR   (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);

extern void   get_BCoord_Exponents (const unsigned int P, const unsigned int d, unsigned int *NExp, unsigned int **NpermsOut,
                                    unsigned int **ExponentsOut);
extern double *basis_TP_Bezier (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double *basis_SI_Bezier (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern void   rst_to_barycentric_SI (const unsigned int Nn, const unsigned int d, const double *rst, double *BCoords);


extern struct S_MATRIX *basis_mat (unsigned int const P, unsigned int const d, unsigned int const Nn,
                                   double const *const rst, basis_tdef basis);

#endif // DPG__bases_h__INCLUDED
