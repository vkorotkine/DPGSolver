// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
/// \file

#include "const_cast.h"

// Standard data types ********************************************************************************************** //

void const_cast_i (const int* dest, const int src)
{
	*(int*) dest = src;
}

void const_cast_ptrdiff (const ptrdiff_t* dest, const ptrdiff_t src)
{
	*(ptrdiff_t*) dest = src;
}

void const_cast_bool (const bool* dest, const bool src)
{
	*(bool*) dest = src;
}

void const_cast_c1 (const char*const* dest, const char*const src)
{
	*(char**) dest = (char*) src;
}

// Custom data types ************************************************************************************************ //
