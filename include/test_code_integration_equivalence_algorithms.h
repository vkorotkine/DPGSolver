// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_equivalence_algorithms_h__INCLUDED
#define DPG__test_code_integration_equivalence_algorithms_h__INCLUDED

#include <stdbool.h>

struct S_equivalence_algs {
	bool         TestTRI;
	char         **argvNew, *PrintName;
	unsigned int P, ML, PG_add, IntOrder_add, IntOrder_mult, NAlgs;
	int          nargc;
};

extern void test_equivalence_algorithms (struct S_equivalence_algs *const data, char const *const TestName);

#endif // DPG__test_code_integration_equivalence_algorithms_h__INCLUDED
