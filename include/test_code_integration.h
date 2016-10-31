// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#ifndef DPG__test_code_integration_h__INCLUDED
#define DPG__test_code_integration_h__INCLUDED

extern void check_convergence_orders (const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
	                                  const unsigned int PMax, unsigned int *pass);
extern void evaluate_mesh_regularity (void);
extern void code_startup             (int nargc, char **argv, const unsigned int Nref, const unsigned int update_argv);
extern void code_cleanup             (void);

#endif // DPG__test_code_integration_h__INCLUDED
