#ifndef DPG__optimize_h__INCLUDED
#define DPG__optimize_h__INCLUDED

struct Simulation;

void optimize(const struct Simulation* sim, const char*const ctrl_file_name);

#endif // DPG__optimize_h__INCLUDED