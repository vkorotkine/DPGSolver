// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)
///	\file

#include "const_cast.h"
#include "Element.h"

// Standard data types ********************************************************************************************** //

void const_cast_ui (const unsigned int* dest, const unsigned int src)
{
	*(unsigned int*) dest = src;
}

void const_cast_st (const size_t* dest, const size_t src)
{
	*(size_t*) dest = src;
}

void const_cast_bool (const bool* dest, const bool src)
{
	*(bool*) dest = src;
}

// Custom data types ************************************************************************************************ //

void const_cast_const_Element (const struct const_Element*const* dest, const struct const_Element*const src)
{
	*(struct const_Element**) dest = (struct const_Element*) src;
}
