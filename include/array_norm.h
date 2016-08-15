// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__array_norm_h__INCLUDED
#define DPG__array_norm_h__INCLUDED

#include "petscmat.h"


extern unsigned int array_norm_ui      (const unsigned int LenA, const unsigned int *A, const char *NormType);
extern double       array_norm_d       (const unsigned int LenA, const double       *A, const char *NormType);
extern unsigned int array_norm_diff_ui (const unsigned int LenA, const unsigned int *A, const unsigned int *B, const char *NormType);
extern double       array_norm_diff_d  (const unsigned int LenA, const double       *A, const double       *B, const char *NormType);

extern double       PetscMatAIJ_norm_diff_d (const unsigned int NRows, Mat A, Mat B, const char *NormType);

#endif // DPG__array_norm_h__INCLUDED
