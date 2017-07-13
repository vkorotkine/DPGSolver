// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Test_h__INCLUDED
#define DPG__Test_h__INCLUDED

#include <stdbool.h>
#include <limits.h>

#define TEST_N_RIEMANN      4
#define TEST_N_LF           2
#define TEST_N_ROE          4
#define TEST_N_BACKPRESSURE 2

#define TEST_N_INVISCID_FLUXES 2
#define TEST_N_VISCOUS_FLUXES  2

struct S_TEST {
	// Counters
	unsigned int Ntest, Npass, Nwarnings,
	             EnteredRiemann[TEST_N_RIEMANN],
	             EnteredLF[TEST_N_LF],
	             EnteredRoe[TEST_N_ROE],
	             EnteredBackPressure[TEST_N_BACKPRESSURE],
	             EnteredInviscidFlux[TEST_N_INVISCID_FLUXES],
	             EnteredViscousFlux[TEST_N_VISCOUS_FLUXES];

	// Initialization
	bool         Active;
	char         *TestCase;
	unsigned int PGlobal, ML;

	// Parameters
	unsigned int PG_add, IntOrder_add, IntOrder_mult;

	// Integration Testing
	bool CheckOffDiagonal;
};
extern struct S_TEST TestDB;

#endif // DPG__Test_h__INCLUDED
