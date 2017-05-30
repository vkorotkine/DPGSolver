// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__cubature_h__INCLUDED
#define DPG__cubature_h__INCLUDED

extern void cubature_TP  (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_TRI (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_TET (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);
extern void cubature_PYR (double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                          const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType);

#include <stdbool.h>

struct S_CUBATURE {
	bool return_w, return_symm;
	char const *NodeType;

	unsigned int Nn, Ns, P, d,
	             *symms;
	double       *rst, *w;
};

extern void cubature_s_TP (struct S_CUBATURE *const CUBDATA);

#endif // DPG__cubature_h__INCLUDED
