// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_code_bases_h__INCLUDED
#define DPG__test_code_bases_h__INCLUDED

extern void get_scaling_basis_TRI (const unsigned int i, const unsigned int j, const double b,
                                   double *con_i, double *con_j, double *con_b);
extern void get_scaling_basis_TET (const unsigned int i, const unsigned int j, const unsigned int k, const double b, const double c,
                                   double *con_i, double *con_j, double *con_k, double *con_b, double *con_c);
extern void get_scaling_basis_PYR (const unsigned int i, const unsigned int j, const unsigned int k, const double c,
                                   double *con_i, double *con_j, double *con_k, double *con_c);
extern void poly2                 (const double *r, const double *s, const double *t, const unsigned int Nn,
                                   double **f, double **f_r, double **f_s, double **f_t);

#endif // DPG__test_code_bases_h__INCLUDED
