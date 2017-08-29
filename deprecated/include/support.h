// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__support_h__INCLUDED
#define DPG__support_h__INCLUDED

#include <stddef.h>
#include <complex.h>

extern void exit_trace (void);
extern void set_to_zero_d     (size_t len, double *A);
extern void set_to_zero_cmplx (size_t len, double complex *A);

#endif // DPG__support_h__INCLUDED
