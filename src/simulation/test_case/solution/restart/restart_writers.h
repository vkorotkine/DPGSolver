/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */

#ifndef DPG__restart_writers_h__INCLUDED
#define DPG__restart_writers_h__INCLUDED
/** \file
 *  \brief Provides the interface to restart file writer containers and functions.
 */

struct Simulation;

/// \brief Container for restart file related information.
struct Restart_Info {
	int ml; ///< 'm'esh 'l'evel for the restart file.
	int p;  ///< 'p'olynomial degree for the restart file.
};

/// \brief Output a restart file for the \ref Simulation.
void output_restart
	(const struct Simulation*const sim,           ///< \ref Simulation.
	 const struct Restart_Info*const restart_info ///< \ref Restart_Info.
	);

#endif // DPG__restart_writers_h__INCLUDED
