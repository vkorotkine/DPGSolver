// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "const_cast.h"

void const_cast_ui (const unsigned int* dest, const unsigned int src)
{
	*(unsigned int*) dest = src;
}

void const_cast_st (const size_t* dest, const size_t src)
{
	*(size_t*) dest = src;
}
