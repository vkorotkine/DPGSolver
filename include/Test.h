// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Test_h__INCLUDED
#define DPG__Test_h__INCLUDED

#include <stdbool.h>

struct S_TEST {
	// Counters
	unsigned int Ntest, Npass, Nwarnings, EnteredRiemann[4], EnteredLF[2], EnteredRoe[4], EnteredBackPressure[2];

	// Initialization
	bool         Active;
	char         *TestCase;
	unsigned int PGlobal, ML;

	// Parameters
	unsigned int PG_add, IntOrder_add, IntOrder_mult;
};
extern struct S_TEST TestDB;

#endif // DPG__Test_h__INCLUDED
