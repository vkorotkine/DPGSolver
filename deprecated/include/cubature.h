// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__cubature_h__INCLUDED
#define DPG__cubature_h__INCLUDED

#include <stdbool.h>

struct S_CUBATURE {
	bool return_w, return_symm;
	char const *NodeType;

	unsigned int Nn, Ns, P, d,
	             *symms;
	double       *rst, *w;
};

typedef void (*cubature_tdef) (struct S_CUBATURE *const CUBDATA);

extern void cubature_TP  (struct S_CUBATURE *const CUBDATA);
extern void cubature_TRI (struct S_CUBATURE *const CUBDATA);
extern void cubature_TET (struct S_CUBATURE *const CUBDATA);
extern void cubature_PYR (struct S_CUBATURE *const CUBDATA);

extern void set_cubdata (struct S_CUBATURE *const CUBDATA, bool const return_w, bool const return_symm,
                         char const *const NodeType, unsigned int const d, unsigned int const P,
                         cubature_tdef cubature);
extern void set_from_cubdata (struct S_CUBATURE const *const CUBDATA, unsigned int *Nn, unsigned int *Ns, double **rst,
                              double **w, unsigned int **symms);

extern struct S_CUBATURE *cub_constructor (bool const return_w, bool const return_symm, char const *const NodeType,
                                           unsigned int const d, unsigned int const P, cubature_tdef cubature);
extern void cub_destructor (struct S_CUBATURE *CUBDATA);

#endif // DPG__cubature_h__INCLUDED
