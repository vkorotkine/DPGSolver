// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "support.h"

#include <stdlib.h>
#include <stdio.h>
#include <execinfo.h>
//#include <signal.h>
//#include <unistd.h>

/*
 *	Purpose:
 *		Provide support functions for the code, generally for debugging purposes: (ToBeModified)
 *			void exit_backtrace(void)
 *
 *	Comments:
 *		This is not currently working, always giving a single stack frame.
 *
 *	Notation:
 *
 *	References:
 *		http://stackoverflow.com/questions/77005/how-to-generate-a-stacktrace-when-my-gcc-c-app-crashes
 *		https://www.gnu.org/software/libc/manual/html_node/Backtraces.html
 *
 */

void exit_trace(void)
{
	unsigned int i;
	void         *array[50];
	size_t       size;
	char         **strings;

	size = backtrace(array,50);
	strings = backtrace_symbols(array,size);

	printf ("Obtained %zd stack frames.\n", size);
	for (i = 0; i < size; i++)
		printf ("%s\n", strings[i]);

	free(strings);
	exit(1);
}
