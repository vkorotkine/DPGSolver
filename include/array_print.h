// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__array_print_h__INCLUDED
#define DPG__array_print_h__INCLUDED

#include <complex.h>


extern void array_print_ui    (const unsigned int m, const unsigned int n, const unsigned int   *A, const char layout);
extern void array_print_i     (const unsigned int m, const unsigned int n, const int            *A, const char layout);
extern void array_print_l     (const unsigned int m, const unsigned int n, const long           *A, const char layout);
extern void array_print_ll    (const unsigned int m, const unsigned int n, const long long      *A, const char layout);
extern void array_print_f     (const unsigned int m, const unsigned int n, const float          *A, const char layout);
extern void array_print_d     (const unsigned int m, const unsigned int n, const double         *A, const char layout);
extern void array_print_ld    (const unsigned int m, const unsigned int n, const long double    *A, const char layout);
extern void array_print_cmplx (const unsigned int m, const unsigned int n, const double complex *A, const char layout);

#endif // DPG__array_print_h__INCLUDED
