// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__setup_periodic_h__INCLUDED
#define DPG__setup_periodic_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>

#include "petscsys.h"
 
#include "Database.h"
#include "Parameters.h"

#include "element_functions.h"
#include "array_find_index.h"
#include "array_print.h"


extern void setup_periodic (void);

#endif // DPG__setup_periodic_h__INCLUDED
