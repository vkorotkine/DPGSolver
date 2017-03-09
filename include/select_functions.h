// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__select_functions_h__INCLUDED
#define DPG__select_functions_h__INCLUDED

typedef void (*cubature_tdef) (double **rst, double **w_vec, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                               const unsigned int return_w, const unsigned int P, const unsigned int d,
                               const char *NodeType);
typedef double *(*basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn, unsigned int *NbfOut,
                               const unsigned int d);
typedef double **(*grad_basis_tdef) (const unsigned int P, const double *rst, const unsigned int Nn,
                                     unsigned int *NbfOut, const unsigned int d);


extern void select_functions            (basis_tdef *basis, grad_basis_tdef *grad_basis, cubature_tdef *cubature,
                                         const unsigned int type);
extern void select_functions_basis      (basis_tdef *basis,           const unsigned int type);
extern void select_functions_grad_basis (grad_basis_tdef *grad_basis, const unsigned int type);
extern void select_functions_cubature   (cubature_tdef *cubature,     const unsigned int type);
extern void select_functions_basis_Bezier (basis_tdef *basis, const unsigned int type);

#endif // DPG__select_functions_h__INCLUDED
