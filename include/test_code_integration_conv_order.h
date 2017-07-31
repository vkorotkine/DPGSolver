// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_conv_order_h__INCLUDED
#define DPG__test_code_integration_conv_order_h__INCLUDED

#include <stdbool.h>

struct S_convorder {
	bool         PrintEnabled, SolveExplicit, SolveImplicit, AdaptiveRefine, omit_root, Compute_L2proj;
	char         **argvNew, *PrintName;
	unsigned int PMin, PMax, MLMin, MLMax, Adapt, PG_add, IntOrder_add, IntOrder_mult;
	int          nargc;
};

extern void test_conv_order (struct S_convorder *const data, char const *const TestName);

#endif // DPG__test_code_integration_conv_order_h__INCLUDED
