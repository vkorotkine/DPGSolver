// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__gmsh_reader_h__INCLUDED
#define DPG__gmsh_reader_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "metis.h"
#include "parmetis.h"
#include "petscsys.h"
#include "mkl.h"
 
#include "Database.h"
#include "Parameters.h"

#include "array_sort.h"
#include "array_print.h"
#include "array_find_index.h"


extern void gmsh_reader               (void);
extern void find_periodic_connections (unsigned int *Pve, unsigned int *pvePointer, const unsigned int VeMax);

#endif // DPG__gmsh_reader_h__INCLUDED
