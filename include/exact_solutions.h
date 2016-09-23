// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__exact_solutions_h__INCLUDED
#define DPG__exact_solutions_h__INCLUDED

extern void compute_exact_solution (const unsigned int Nn, double *XYZ, double *UEx, double *sEx, const unsigned int solved);
extern void compute_source         (const unsigned int Nn, double *XYZ, double *source);

#endif // DPG__exact_solutions_h__INCLUDED
