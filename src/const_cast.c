// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "const_cast.h"

void const_cast_ui_1 (const unsigned int*const* dest, unsigned int* src)
{
	*(unsigned int**)& dest = src;
}

void const_cast_ui_2 (const unsigned int*const*const* dest, unsigned int** src)
{
	*(unsigned int***)& dest = src;
}

void const_cast_st_0 (const size_t* dest, const size_t src)
{
	*(size_t*)& dest = src;
}
