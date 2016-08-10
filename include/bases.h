// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__bases_h__INCLUDED
#define DPG__bases_h__INCLUDED

extern double *basis_TP        (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double *basis_SI        (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double *basis_PYR       (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_TP  (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_SI  (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern double **grad_basis_PYR (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut, const unsigned int d);
extern void   rst_to_abc_SI    (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);
extern void   rst_to_abc_PYR   (const unsigned int Nn, const unsigned int d, const double *rst, double *a, double *b, double *c);

#endif // DPG__bases_h__INCLUDED
