// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_h__INCLUDED
#define DPG__test_code_integration_h__INCLUDED

#include <stdbool.h>

extern void code_startup             (int const nargc, char const *const *const argv, unsigned int const Nref,
                                      unsigned int const update_argv);
extern void code_cleanup             (void);
extern void code_startup_mod_prmtrs  (int const nargc, char const *const *const argv, unsigned int const Nref,
                                      unsigned int const update_argv, unsigned int const phase);
extern void code_startup_mod_ctrl    (int const nargc, char const *const *const argv, unsigned int const Nref,
                                      unsigned int const update_argv, unsigned int const phase);
extern void evaluate_mesh_regularity (double *mesh_quality);
extern void check_convergence_orders (const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                                      const unsigned int PMax, unsigned int *pass, const bool PrintEnabled);
extern void check_mesh_regularity    (const double *mesh_quality, const unsigned int NML, unsigned int *pass,
                                      const bool PrintEnabled);
extern void set_PrintName            (char *name_type, char *PrintName, bool *omit_root);
extern void compute_dof              (void);

#endif // DPG__test_code_integration_h__INCLUDED
