// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__cubature_h__INCLUDED
#define DPG__cubature_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "petscsys.h"
#include "mkl.h"
 
#include "Parameters.h"

#include "array_norm.h"
#include "array_sort.h"
#include "matrix_functions.h"


extern void cubature_TP  (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_TRI (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_TET (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_PYR (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);

#endif // DPG__cubature_h__INCLUDED
