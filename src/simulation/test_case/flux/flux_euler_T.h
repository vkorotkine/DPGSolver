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
/** \file
 *  \brief Provides templated functions relating to euler fluxes.
 */

struct Flux_Input_T;
struct mutable_Flux_T;

/** \brief Version of \ref compute_Flux_fptr computing the fluxes for the Euler equations.
 *
 *  Euler fluxes (eq, dim):
 *  \f[
 *  f =
 *  \begin{bmatrix}
 *  	\rho u & \rho v & \rho w \\
 *  	\rho u^2 + p & \rho uv & \rho uw \\
 *  	\rho vu & \rho v^2 + p & \rho vw \\
 *  	\rho wu & \rho wv & \rho w^2 + p \\
 *  	(E+p)u & (E+p)v & (E+p)w \\
 *  \end{bmatrix}
 *  \f]
 */
void compute_Flux_T_euler
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);

/** \brief Version of \ref compute_Flux_fptr computing the fluxes and Jacobians for the Euler equations.
 *
 *  Euler flux Jacobians (eq, var, dim):
 *  \f[
 *  \frac{df}{ds} =
 *  \begin{bmatrix}
 *  \begin{bmatrix}
 *  	 0                          &  1               &  0               &  0               & 0        \\
 *  	-u^2+\frac{\gamma-1}{2}V^2  & -(\gamma-3)u     & -(\gamma-1)v     & -(\gamma-1)w     & \gamma-1 \\
 *  	-uv                         &  v               &  u               &  0               & 0        \\
 *  	-uw                         &  w               &  0               &  u               & 0        \\
 *  	 u(\frac{\gamma-1}{2}V^2-H) &  H-(\gamma-1)u^2 & -(\gamma-1)uv    & -(\gamma-1)uw    & \gamma u \\
 *  \end{bmatrix},
 *  \begin{bmatrix}
 *  	 0                          &  0               &  1               &  0               & 0        \\
 *  	-uv                         &  v               &  u               &  0               & 0        \\
 *  	-v^2+\frac{\gamma-1}{2}V^2  & -(\gamma-1)u     & -(\gamma-3)v     & -(\gamma-1)w     & \gamma-1 \\
 *  	-vw                         &  0               &  w               &  v               & 0        \\
 *  	 v(\frac{\gamma-1}{2}V^2-H) & -(\gamma-1)uv    &  H-(\gamma-1)v^2 & -(\gamma-1)vw    & \gamma v \\
 *  \end{bmatrix},
 *  \begin{bmatrix}
 *  	 0                          &  0               &  0               &  1               & 0        \\
 *  	-uw                         &  w               &  0               &  u               & 0        \\
 *  	-vw                         &  0               &  w               &  v               & 0        \\
 *  	-w^2+\frac{\gamma-1}{2}V^2  & -(\gamma-1)u     & -(\gamma-1)v     & -(\gamma-3)w     & \gamma-1 \\
 *  	 w(\frac{\gamma-1}{2}V^2-H) & -(\gamma-1)uw    & -(\gamma-1)vw    &  H-(\gamma-1)w^2 & \gamma w \\
 *  \end{bmatrix},
 *  \end{bmatrix}
 *  \f]
 */
void compute_Flux_T_euler_jacobian
	(const struct Flux_Input_T* flux_i, ///< See brief.
	 struct mutable_Flux_T* flux        ///< See brief.
	);
