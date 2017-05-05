// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_boundary_conditions_h__INCLUDED
#define DPG__test_code_boundary_conditions_h__INCLUDED

#include <stdbool.h>

extern void set_memory_test_boundary_conditions     (char const operation);
extern void set_BTypes                              (unsigned int *NBTypesOut, char ***BTypeOut);
extern void set_parameters_test_boundary_conditions (char const *const BType, unsigned int const d);
extern void reset_entered_test_boundary_conditions  (char const *const BType);
extern void update_values_BackPressure              (unsigned int const Nn, unsigned int const Nel, double *const W,
                                                     double *const nL, unsigned int const d);
extern void check_entered_test_boundary_conditions  (bool *CheckedAll, char const *const BType);

#endif // DPG__test_code_boundary_conditions_h__INCLUDED
