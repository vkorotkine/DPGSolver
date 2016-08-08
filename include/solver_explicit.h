// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__solver_explicit_h__INCLUDED
#define DPG__solver_explicit_h__INCLUDED

#include <stdlib.h>
#include <stdio.h>
#include <math.h> // ToBeModified
#include <string.h>
 
#include "Database.h"
#include "Parameters.h"
#include "Macros.h"

#include "adaptation.h"
#include "update_VOLUMEs.h"
#include "explicit_VOLUME_info.h"
#include "explicit_FACET_info.h"
#include "finalize_RHS.h"
#include "output_to_paraview.h"


extern void solver_explicit (void);

#endif // DPG__solver_explicit_h__INCLUDED
