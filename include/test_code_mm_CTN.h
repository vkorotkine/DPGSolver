// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_mm_CTN_h__INCLUDED
#define DPG__test_code_mm_CTN_h__INCLUDED

extern void mm_CTN_mv_fully_unrolled_d (const int m, const int n, const int k, double *A, double *B, double *C);
extern void mm_CTN_mvBLAS_d            (const int m, const int n, const int k, double *A, double *B, double *C);
extern void mm_CTN_fully_unrolled_mv_d (const int m, const int n, const int k, double *A, double *B, double *C);
extern void mm_CTN_mmBLAS_d            (const int m, const int n, const int k, double *A, double *B, double *C);

#endif // DPG__test_code_mm_CTN_h__INCLUDED
