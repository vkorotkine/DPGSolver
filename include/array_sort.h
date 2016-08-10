// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__array_sort_h__INCLUDED
#define DPG__array_sort_h__INCLUDED

extern void array_sort_ui (unsigned int NRows, unsigned int NCols, unsigned int *A, unsigned int *Indices, const char ordering, const char trans);
extern void array_sort_i  (unsigned int NRows, unsigned int NCols, int          *A, unsigned int *Indices, const char ordering, const char trans);
extern void array_sort_d  (unsigned int NRows, unsigned int NCols, double       *A, unsigned int *Indices, const char ordering, const char trans);

#endif // DPG__array_sort_h__INCLUDED
