#include <stdio.h>

#include "test.h"

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

void test_print(const unsigned int pass)
{
	TestDB.Ntest++;

	if (pass)
		printf("Pass\n");
	else
		printf("Fail --- Fail --- Fail\n");
}
