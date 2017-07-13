// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_support.h"

#include <stdio.h>
#include <stdlib.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "Test.h"

void test_print2(unsigned int const pass, char const *const string)
{
	/*
	 *	Purpose:
	 *		Print pass/failure message for testing functions.
	 *
	 *	Comments:
	 *
	 *	Notation:
	 *
	 *	References:
	 */

	TestDB.Ntest++;

	printf("%-90s",string);
	if (pass) {
		TestDB.Npass++;
		printf("Pass\n");
	} else {
		printf("Fail --- Fail --- Fail\n");
	}
}

void test_print_warning(char const *const string)
{
	TestDB.Nwarnings++;

	printf("\n\n********************************************************************************************\n");
	printf("Warning: %s.\n",string);
	printf("********************************************************************************************\n\n\n");
}

void set_memory_test_jacobians(char const operation)
{
	if (operation == 'A') { // (A)llocate
		DB.TestCase      = calloc(STRLEN_MAX , sizeof *(DB.TestCase));      // free
		DB.PDE           = calloc(STRLEN_MAX , sizeof *(DB.PDE));           // free
		DB.PDESpecifier  = calloc(STRLEN_MAX , sizeof *(DB.PDESpecifier));  // free
		DB.Geometry      = calloc(STRLEN_MAX , sizeof *(DB.Geometry));      // free
		DB.GeomSpecifier = calloc(STRLEN_MAX , sizeof *(DB.GeomSpecifier)); // free
		DB.MeshFile      = calloc(STRLEN_MAX , sizeof *(DB.MeshFile));      // free
		DB.MeshType      = calloc(STRLEN_MAX , sizeof *(DB.MeshType));      // free
	} else if (operation == 'F') { // (F)ree
		free(DB.TestCase);
		free(DB.PDE);
		free(DB.PDESpecifier);
		free(DB.Geometry);
		free(DB.GeomSpecifier);
		free(DB.MeshFile);
		free(DB.MeshType);
	} else {
		EXIT_UNSUPPORTED;
	}
}
