// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_integration_h__INCLUDED
#define DPG__test_code_integration_h__INCLUDED
/// \file

#include <stdbool.h>

void code_startup
	(const int nargc,                 ///< Standard.
	 const char*const*const argv,     ///< Standard.
	 const unsigned int n_ref,        ///< (n)umber of (ref)inement steps to perform.
	 const unsigned int modify_params /**< Flag for how parameters have been modified.\n
	                                   *   Options:
	                                   *   - 0: No modifications;
	                                   *   - 1: Modified order;
	                                   *   - 2: Modified order and mesh level;
	                                   */
	);

extern void code_cleanup             (void);
extern void code_startup_mod_prmtrs  (const int nargc, char const *const *const argv, const unsigned int Nref,
                                      const unsigned int update_argv, const unsigned int phase);
extern void code_startup_mod_ctrl    (const int nargc, char const *const *const argv, const unsigned int Nref,
                                      const unsigned int update_argv, const unsigned int phase);
extern void evaluate_mesh_regularity (double *mesh_quality);
extern void check_convergence_orders (const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                                      const unsigned int PMax, unsigned int *pass, const bool PrintEnabled);
extern void check_mesh_regularity    (const double *mesh_quality, const unsigned int NML, unsigned int *pass,
                                      const bool PrintEnabled);
extern void set_PrintName            (char *name_type, char *PrintName, bool *omit_root);
extern void compute_dof              (void);

#endif // DPG__test_code_integration_h__INCLUDED
