// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__variable_functions_h__INCLUDED
#define DPG__variable_functions_h__INCLUDED

extern void convert_variables (double *VarIn, double *VarOut, const unsigned int dIn, const unsigned int dOut,
                               const unsigned int Nn, const unsigned int Nel, const char TypeIn, const char TypeOut);

#endif // DPG__variable_functions_h__INCLUDED
