// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "exact_solutions.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions related to exact solutions.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */
static double generate_Euler_source(unsigned int eq, double *RUVP) ;

static double get_boundary_value_Advection (double const x, double const y, double const z);

void compute_exact_solution(const unsigned int Nn, const double *XYZ, double *UEx, const unsigned int solved)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i;
	double       *rhoEx, *uEx, *vEx, *wEx, *pEx;
	const double *X, *Y, *Z;

	rhoEx = &UEx[Nn*0];
	uEx   = &UEx[Nn*1];
	vEx   = &UEx[Nn*2];
	wEx   = &UEx[Nn*3];
	pEx   = &UEx[Nn*4];

	X = &XYZ[0*Nn];
	Y = &XYZ[1*Nn];
	Z = &XYZ[(d-1)*Nn];

	if (strstr(TestCase,"PeriodicVortex")) {
		// Initialize DB Parameters
		double       pInf           = DB.pInf,
					 TInf           = DB.TInf,
					 uInf           = DB.uInf,
					 vInf           = DB.vInf,
					 wInf           = DB.wInf,
					 Rg             = DB.Rg,
					 Cscale         = DB.Cscale,
					 PeriodL        = DB.PeriodL,
					 PeriodFraction = DB.PeriodFraction,
					 Xc0            = DB.Xc,
					 Yc             = DB.Yc,
					 Rc             = DB.Rc;

		// Standard datatypes
		double DistTraveled, Xc, rhoInf, r2, C;

		rhoInf = pInf/(Rg*TInf);
		C      = Cscale;

		if (solved) {
			DistTraveled = PeriodL*PeriodFraction;
			Xc = Xc0 + DistTraveled;
			while (Xc > 0.5*PeriodL)
				Xc -= PeriodL;
		} else {
			Xc = Xc0;
		}

		for (i = 0; i < Nn; i++) {
			r2 = (pow(X[i]-Xc,2.0)+pow(Y[i]-Yc,2.0))/(Rc*Rc);

			uEx[i]   = uInf - C*(Y[i]-Yc)/(Rc*Rc)*exp(-0.5*r2);
			vEx[i]   = vInf + C*(X[i]-Xc)/(Rc*Rc)*exp(-0.5*r2);
			wEx[i]   = wInf;
			pEx[i]   = pInf - rhoInf*(C*C)/(2*Rc*Rc)*exp(-r2);
			rhoEx[i] = rhoInf;
		}
	} else if (strstr(TestCase,"SupersonicVortex")) {
		// Initialize DB Parameters
		double rIn   = DB.rIn,
		       MIn   = DB.MIn,
		       rhoIn = DB.rhoIn,
		       VIn   = DB.VIn;

		// Standard datatypes
		double r, t, Vt;

		for (i = 0; i < Nn; i++) {
			r = sqrt(X[i]*X[i]+Y[i]*Y[i]);
			t = atan2(Y[i],X[i]);

			rhoEx[i] = rhoIn*pow(1.0+0.5*GM1*MIn*MIn*(1.0-pow(rIn/r,2.0)),1.0/GM1);
			pEx[i]   = pow(rhoEx[i],GAMMA)/GAMMA;

			Vt = -VIn/r;
			uEx[i] = -sin(t)*Vt;
			vEx[i] =  cos(t)*Vt;
			wEx[i] =  0.0;
		}
	} else if (strstr(TestCase,"PlaneCouette")) {
		double const uIn = DB.uIn,
		             TIn = DB.TIn,
		             pIn = DB.pIn,
		             Pr  = DB.Pr,
		             Rg  = DB.Rg;

		double const b = uIn*uIn/(DB.Cp*TIn);

		if (DB.Const_mu) {
			for (size_t n = 0; n < Nn; n++) {
				double const eta = 0.5*(Y[n]+1.0),
				             T   = -0.5*Pr*b*pow(1.0-eta,2.0)+1.0+0.5*Pr*b;

				pEx[n]   = pIn;
				rhoEx[n] = pEx[n]/(Rg*T);
				uEx[n]   = 1.0-eta;
				vEx[n]   = 0.0;
				wEx[n]   = 0.0;
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestCase,"TaylorCouette")) {
		// Note: This exact solution is valid only before the Taylor-Couette instability develops and is only accurate
		//       for velocity and temperature (except for r = rIn where all components are exact).
		if (d != 2)
			EXIT_UNSUPPORTED;

		double rIn   = DB.rIn,
		       rOut  = DB.rOut,
		       omega = DB.omega,
		       TIn   = DB.TIn,
		       pIn   = DB.pIn,
		       mu    = DB.mu,
		       kappa = DB.kappa,
		       Rg    = DB.Rg;

		double C = omega/(1.0/(rIn*rIn)-1.0/(rOut*rOut));

		for (size_t n = 0; n < Nn; n++) {
			double r, t, Vt, T;
			r = sqrt(X[n]*X[n]+Y[n]*Y[n]);
			t = atan2(Y[n],X[n]);

			Vt = C*(1.0/r-r/(rOut*rOut));
			T  = TIn - 2.0*C*C/(rOut*rOut)*mu/kappa*log(r/rIn) - C*C*mu/kappa*(1.0/(r*r)-1.0/(rIn*rIn));

			// Illingworth(1950), p.8 notes that the pressure is nearly uniform => set pEx ~= pIn and compute rhoEx
			// using the ideal gas law.
			pEx[n]   = pIn;
			rhoEx[n] = pEx[n]/(Rg*T);
			uEx[n]   = -sin(t)*Vt;
			vEx[n]   =  cos(t)*Vt;
			wEx[n]   = 0.0;
		}
	} else if (strstr(TestCase,"Advection")) {
		if (DB.SourcePresent) {
			EXIT_UNSUPPORTED; // Manufactured solution
		} else {
			if (strstr(DB.PDESpecifier,"Unsteady")) {
				EXIT_UNSUPPORTED;
			} else if (strstr(DB.PDESpecifier,"Steady")) {
				// Extrapolate from upwind boundary
				for (size_t n = 0; n < Nn; n++)
					UEx[n] = get_boundary_value_Advection(X[n],Y[n],Z[n]);
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	} else if (strstr(TestCase,"Poisson")) {
		double Poisson_scale = DB.Poisson_scale;

		if (fabs(Poisson_scale) < EPS)
			printf("Error: Make sure to set Poisson_scale.\n"), EXIT_MSG;

		for (i = 0; i < Nn; i++) {
			if (d == 1)
				UEx[i] = sin(2.0*PI*X[i]);
			else if (d == 2)
//				UEx[i] = sin(PI*X[i])*sin(PI*Y[i]);
				UEx[i] = cos(Poisson_scale*PI*X[i])*cos(Poisson_scale*PI*Y[i]);
//				UEx[i] = X[i]*Y[i]*sin(PI*X[i])*sin(PI*Y[i]);
			else if (d == 3)
				UEx[i] = sin(PI*X[i])*sin(PI*Y[i])*sin(PI*Z[i]);
			else
				EXIT_UNSUPPORTED;
		}
	} else if (strstr(TestCase,"EllipticPipe")) {

		double rho_0 = DB.rho_store[0],
		       rho_1 = DB.rho_store[1],
		       rho_2 = DB.rho_store[2],
		       rho_3 = DB.rho_store[3];
		       //rho_2 = DB.rho_store[4],
		       //rho_xy = DB.rho_store[5],
		       //rho_3 = DB.rho_store[6];

		double p_0 = DB.p_store[0],
		       p_1 = DB.p_store[1],
		       p_2 = DB.p_store[2],
		       p_3 = DB.p_store[3];
		       //p_2 = DB.p_store[4],
		       //p_xy = DB.p_store[5],
		       //p_3 = DB.p_store[6];

		double u_0 = DB.u_store[0],
		       u_1 = DB.u_store[1],
		       u_2 = DB.u_store[2],
		       u_3 = DB.u_store[3];
		       //u_2 = DB.u_store[4],
		       //u_xy = DB.u_store[5],
		       //u_3 = DB.u_store[6];

		double v_0 = DB.v_store[0],
		       v_1 = DB.v_store[1],
		       v_2 = DB.v_store[2],
		       v_3 = DB.v_store[3];
		       //v_2 = DB.v_store[4],
		       //v_xy = DB.v_store[5],
		       //v_3 = DB.v_store[6];

		/*double f_0 = DB.f_store[0],
		       f_1 = DB.f_store[1],
		       f_2 = DB.f_store[2],
		       f_3 = DB.f_store[3],
		       f_4 = DB.f_store[4],
		       f_5 = DB.f_store[5],
		       f_6 = DB.f_store[6];*/

		/*double a = DB.geo_store[0],
			   b = DB.geo_store[1],
			   c = DB.geo_store[2];*/

		//MMS1
		/*for (size_t n = 0; n < Nn; n++) {
			rhoEx[n] = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			pEx[n] = p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			uEx[n] = u_0 + u_x*sin(u_1*X[n])+u_y*sin(u_2*Y[n])+u_xy*sin(u_3*X[n]*Y[n]);
			vEx[n] = v_0 + v_x*sin(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*sin(v_3*X[n]*Y[n]);
			wEx[n] = 0.0;
			}*/

		//MMS2
		/*for (size_t n = 0; n < Nn; n++) {
			rhoEx[n] = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n])+rho_xy*cos(rho_3*X[n]*Y[n]);
			pEx[n] = p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			uEx[n] = u_0 + u_x*sin(u_1*X[n])+u_y*cos(u_2*Y[n])+u_xy*cos(u_3*X[n]*Y[n]);
			vEx[n] = v_0 + v_x*cos(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*cos(v_3*X[n]*Y[n]);
			wEx[n]   = 0.0;
		}*/

		//MMS3
		/*for (size_t n = 0; n < Nn; n++) {
			double const xi = (2/c)*X[n]-1, //The domain transformations
						 eta = (2*a*Y[n])/(b*sqrt(a*a+X[n]*X[n]))-1;

			double const u_a = a*sqrt(a*a+0.25*c*c*(xi+1)*(xi+1)),
						 u_b = u_1*sin(u_2*eta)+u_3*cos(u_4*eta),
						 u_c = f_0 + f_1*sin(f_2*(xi+1))+f_3*sin(f_4*(eta+1))+f_5*sin(f_6*(xi+1)*(eta+1)),
						 v_a = 0.5*b*c*(xi+1),
						 v_b = v_1*sin(v_2*eta)+v_3*cos(v_4*eta),
						 v_c = f_0 + f_1*sin(f_2*(xi+1))+f_3*sin(f_4*(eta+1))+f_5*sin(f_6*(xi+1)*(eta+1));

			rhoEx[n] = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			pEx[n] = p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			uEx[n] = u_a*u_b*u_c;
			vEx[n] = v_a*v_b*v_c;
			wEx[n] = 0.0;
		}*/

		//MMS4
		for (size_t n = 0; n < Nn; n++) {
			rhoEx[n] = rho_0 + rho_1*cos(rho_2*X[n])*cos(rho_3*Y[n]);
			pEx[n] = p_0 + p_1*cos(p_2*X[n])*cos(p_3*Y[n]);
			uEx[n] = u_0 + u_1*cos(u_2*X[n])*cos(u_3*Y[n]);
			vEx[n] = v_0 + v_1*cos(v_2*X[n])*cos(v_3*Y[n]);
			wEx[n] = 0.0;
		}

	} else if (strstr(TestCase,"ParabolicPipe")) {
		if (d != 2)
			EXIT_UNSUPPORTED;

		double rho_0 = DB.rho_store[0],
		       rho_x = DB.rho_store[1],
		       rho_1 = DB.rho_store[2],
		       rho_y = DB.rho_store[3],
		       rho_2 = DB.rho_store[4],
		       rho_xy = DB.rho_store[5],
		       rho_3 = DB.rho_store[6];

		double p_0 = DB.p_store[0],
		       p_x = DB.p_store[1],
		       p_1 = DB.p_store[2],
		       p_y = DB.p_store[3],
		       p_2 = DB.p_store[4],
		       p_xy = DB.p_store[5],
		       p_3 = DB.p_store[6];

		double u_0 = DB.u_store[0],
		       u_x = DB.u_store[1],
		       u_1 = DB.u_store[2],
		       u_y = DB.u_store[3],
		       u_2 = DB.u_store[4],
		       u_xy = DB.u_store[5],
		       u_3 = DB.u_store[6];

	 	double v_0 = DB.v_store[0],
		       v_x = DB.v_store[1],
		       v_1 = DB.v_store[2],
		       v_y = DB.v_store[3],
		       v_2 = DB.v_store[4],
      	       v_xy = DB.v_store[5],
		       v_3 = DB.v_store[6];

		/*double a_1 = DB.geo_store[0],
			   a_2 = DB.geo_store[1],
			   b = DB.geo_store[2];*/

		/*for (size_t n = 0; n < Nn; n++) {
			rhoEx[n] = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n])+rho_xy*cos(rho_3*X[n]*Y[n]);
			pEx[n] = p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			uEx[n] = u_0 + u_x*sin(u_1*X[n])+u_y*cos(u_2*Y[n])+u_xy*cos(u_3*X[n]*Y[n]);
			vEx[n] = v_0 + v_x*cos(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*cos(v_3*X[n]*Y[n]);
			wEx[n]   = 0.0;
		}*/

		for (size_t n = 0; n < Nn; n++) {
			rhoEx[n] = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			pEx[n] = p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			uEx[n] = u_0 + u_x*sin(u_1*X[n])+u_y*sin(u_2*Y[n])+u_xy*sin(u_3*X[n]*Y[n]);
			vEx[n] = v_0 + v_x*sin(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*sin(v_3*X[n]*Y[n]);
			wEx[n] = 0.0;
		}

		/*for (size_t n = 0; n < Nn; n++) {

			double const xi = 2*sqrt(b/a_1)*X[n]-1,
						 eta = (2*Y[n]+2*b*X[n]*X[n]-a_1-a_2)/(a_2-a_1);

			rhoEx[n] = rho_0 + rho_1*sin(rho_2*(xi+rho_3))+rho_4*sin(rho_5*(eta+rho_6));
			pEx[n] = p_0 + p_1*sin(p_2*(xi+p_3))+p_4*sin(p_5*(eta+p_6));
			uEx[n] = u_0 + u_1*sin(u_2*(xi+u_3))+u_4*sin(u_5*(eta+u_6));
			vEx[n] = v_0 + v_1*sin(v_2*(xi+v_3))+v_4*sin(v_5*(eta+v_6));
			wEx[n] = 0.0;
		}*/



		} else if (strstr(TestCase,"SinusoidalPipe")) {
                    double rho_0 = DB.rho_store[0],
                           rho_x = DB.rho_store[1],
                           rho_1 = DB.rho_store[2],
                           rho_y = DB.rho_store[3],
                           rho_2 = DB.rho_store[4];

                    double p_0 = DB.p_store[0],
                           p_x = DB.p_store[1],
                           p_1 = DB.p_store[2],
                           p_y = DB.p_store[3],
                           p_2 = DB.p_store[4];

                    double w_0 = DB.w_store[0],
                           w_x = DB.w_store[1],
                           w_1 = DB.w_store[2],
                           w_y = DB.w_store[3],
                           w_2 = DB.w_store[4];

                    double a = 1, b = 2;

                    if (d == 2 ) {
                         for (i = 0; i < Nn; i++) {
                      rhoEx[i] = rho_0 + rho_x*sin(rho_1*X[i])+rho_y*cos(rho_2*Y[i]);
                      pEx[i] = p_0 + p_x*cos(p_1*X[i])+p_y*sin(p_2*Y[i]);
                      uEx[i] = w_0 + w_x*sin(w_1*X[i])+w_y*cos(w_2*Y[i]);
                      vEx[i] = (-a*b*sin(b*X[i]))*(w_0 + w_x*sin(w_1*X[i])+w_y*cos(w_2*Y[i]));
                      wEx[i] = 0.0;
                         }
                    } else
                            EXIT_UNSUPPORTED;

        } else {
		EXIT_UNSUPPORTED;
	}
}

void compute_exact_gradient(const unsigned int Nn, const double *XYZ, double *QEx)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i;
	const double *X, *Y, *Z;

	X = &XYZ[0*Nn];
	Y = &XYZ[1*Nn];
	Z = &XYZ[(d-1)*Nn];

	if (strstr(TestCase,"Poisson")) {
		double Poisson_scale = DB.Poisson_scale;
		for (i = 0; i < Nn; i++) {
			if (d == 1) {
				QEx[Nn*0+i] = 2.0*PI*cos(2.0*PI*X[i]);
			} else if (d == 2) {
//				QEx[Nn*0+i] = PI*cos(PI*X[i])*sin(PI*Y[i]);
//				QEx[Nn*1+i] = PI*sin(PI*X[i])*cos(PI*Y[i]);
				QEx[Nn*0+i] = -Poisson_scale*PI*sin(Poisson_scale*PI*X[i])*cos(Poisson_scale*PI*Y[i]);
				QEx[Nn*1+i] = -Poisson_scale*PI*cos(Poisson_scale*PI*X[i])*sin(Poisson_scale*PI*Y[i]);
//				QEx[Nn*0+i] = Y[i]*sin(PI*Y[i])*(sin(PI*X[i])+X[i]*PI*cos(PI*X[i]));
//				QEx[Nn*1+i] = X[i]*sin(PI*X[i])*(sin(PI*Y[i])+Y[i]*PI*cos(PI*Y[i]));
			} else if (d == 3) {
				QEx[Nn*0+i] = PI*cos(PI*X[i])*sin(PI*Y[i])*sin(PI*Z[i]);
				QEx[Nn*1+i] = PI*sin(PI*X[i])*cos(PI*Y[i])*sin(PI*Z[i]);
				QEx[Nn*2+i] = PI*sin(PI*X[i])*sin(PI*Y[i])*cos(PI*Z[i]);
			} else {
				EXIT_UNSUPPORTED;
			}
		}
	} else if (strstr(TestCase,"TaylorCouette")) {
		// Return gradients of velocity components and temperature.
		double const rIn   = DB.rIn,
		             rOut  = DB.rOut,
		             omega = DB.omega,
		             mu    = DB.mu,
		             kappa = DB.kappa;

		double const C = omega/(1.0/(rIn*rIn)-1.0/(rOut*rOut));
		if (d == 2) {
			for (size_t n = 0; n < Nn; n++) {
				double const x = X[n],
				             y = Y[n],
				             r = sqrt(x*x+y*y),
				             t = atan2(y,x);

				double const drdX[DMAX] = {x/r, y/r, 0.0 },
				             dtdX[DMAX] = {-y/(r*r), x/(r*r), 0.0 };

				double const Vt = C*(1.0/r-r/(rOut*rOut));

				double const dVtdr = C*(-1.0/(r*r)-1.0/(rOut*rOut)),
				             dTdr  = -C*C*mu/kappa*(2.0/r*(1.0/(rOut*rOut)-1.0/(r*r)));

				size_t dim = 0, var = 0;
				QEx[Nn*(dim*3+var++)+n] = -(sin(t)*(dVtdr*drdX[dim]) + cos(t)*dtdX[dim]*Vt);
				QEx[Nn*(dim*3+var++)+n] =  (cos(t)*(dVtdr*drdX[dim]) - sin(t)*dtdX[dim]*Vt);
				QEx[Nn*(dim*3+var++)+n] =  (dTdr*drdX[dim]);

				dim++; var = 0;
				QEx[Nn*(dim*3+var++)+n] = -(sin(t)*(dVtdr*drdX[dim]) + cos(t)*dtdX[dim]*Vt);
				QEx[Nn*(dim*3+var++)+n] =  (cos(t)*(dVtdr*drdX[dim]) - sin(t)*dtdX[dim]*Vt);
				QEx[Nn*(dim*3+var++)+n] =  (dTdr*drdX[dim]);
			}
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		EXIT_UNSUPPORTED;
	}
}

void compute_source(const unsigned int Nn, const double *XYZ, double *source)
{
	/*
	 *	Purpose:
	 *		Computes source terms.
	 */

	// Initialize DB Parameters
	char         *TestCase = DB.TestCase;
	unsigned int d         = DB.d,
	             Neq       = DB.Neq;

	// Standard datatypes
	unsigned int n, eq;

	const double *const X = &XYZ[Nn*0],
	             *const Y = &XYZ[Nn*1],
	             *const Z = &XYZ[Nn*(d-1)];


	if (strstr(TestCase,"Poisson")) {
		double Poisson_scale = DB.Poisson_scale;

		for (eq = 0; eq < Neq; eq++) {
			for (n = 0; n < Nn; n++) {
				if (d == 1)
					source[eq*Nn+n] = -pow(2.0*PI,2.0)*sin(2.0*PI*X[n]);
				else if (d == 2)
//					source[eq*Nn+n] = -2.0*PI*PI*sin(PI*X[n])*sin(PI*Y[n]);
					source[eq*Nn+n] = -2.0*pow(Poisson_scale*PI,2.0)*cos(Poisson_scale*PI*X[n])*cos(Poisson_scale*PI*Y[n]);
//					source[eq*Nn+n] = PI*(Y[n]*sin(PI*Y[n])*(2*cos(PI*X[n])-PI*X[n]*sin(PI*X[n]))
//					                     +X[n]*sin(PI*X[n])*(2*cos(PI*Y[n])-PI*Y[n]*sin(PI*Y[n])));
				else if (d == 3)
					source[eq*Nn+n] = -3.0*PI*PI*sin(PI*X[n])*sin(PI*Y[n])*sin(PI*Z[n]);
				else
					EXIT_UNSUPPORTED;
				source[eq*Nn+n] *= -1.0; // Conflict between source term treatment for Euler and Poisson (ToBeDeleted)
			}
		}
	} else if (strstr(TestCase,"EllipticPipe")) {

		if (d != 2)
			EXIT_UNSUPPORTED;

		double rho_0 = DB.rho_store[0],
		       rho_1 = DB.rho_store[1],
		       rho_2 = DB.rho_store[2],
		       rho_3 = DB.rho_store[3];
		       //rho_2 = DB.rho_store[4],
		       //rho_xy = DB.rho_store[5],
		       //rho_3 = DB.rho_store[6];

		double p_0 = DB.p_store[0],
		       p_1 = DB.p_store[1],
		       p_2 = DB.p_store[2],
		       p_3 = DB.p_store[3];
		       //p_2 = DB.p_store[4],
		       //p_xy = DB.p_store[5],
		       //p_3 = DB.p_store[6];

		double u_0 = DB.u_store[0],
		       u_1 = DB.u_store[1],
		       u_2 = DB.u_store[2],
		       u_3 = DB.u_store[3];
		       //u_2 = DB.u_store[4],
		       //u_xy = DB.u_store[5],
		       //u_3 = DB.u_store[6];

		double v_0 = DB.v_store[0],
		       v_1 = DB.v_store[1],
		       v_2 = DB.v_store[2],
		       v_3 = DB.v_store[3];
		       //v_2 = DB.v_store[4],
		       //v_xy = DB.v_store[5],
		       //v_3 = DB.v_store[6];

		/*double f_0 = DB.f_store[0],
		       f_1 = DB.f_store[1],
		       f_2 = DB.f_store[2],
		       f_3 = DB.f_store[3],
		       f_4 = DB.f_store[4],
		       f_5 = DB.f_store[5],
		       f_6 = DB.f_store[6];*/

		/*double a = DB.geo_store[0],
			   b = DB.geo_store[1],
			   c = DB.geo_store[2];*/

		double RUVP[12];

		//MMS1
		/*for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			RUVP[1]  =  rho_x*rho_1*cos(rho_1*X[n])+rho_xy*rho_3*Y[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[2]  =  rho_y*rho_2*cos(rho_2*Y[n])+rho_xy*rho_3*X[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[3]  =  u_0 + u_x*sin(u_1*X[n])+u_y*sin(u_2*Y[n])+u_xy*sin(u_3*X[n]*Y[n]);
			RUVP[4]  =  u_x*u_1*cos(u_1*X[n])+u_xy*u_3*Y[n]*cos(u_3*X[n]*Y[n]);
			RUVP[5]  =  u_y*u_2*cos(u_2*Y[n])+u_xy*u_3*X[n]*cos(u_3*X[n]*Y[n]);
			RUVP[6]  =  v_0 + v_x*sin(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*sin(v_3*X[n]*Y[n]);
			RUVP[7]  =  v_x*v_1*cos(v_1*X[n])+v_xy*v_3*Y[n]*cos(v_3*X[n]*Y[n]);
			RUVP[8]  =  v_y*v_2*cos(v_2*Y[n])+v_xy*v_3*X[n]*cos(v_3*X[n]*Y[n]);
			RUVP[9]  =  p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			RUVP[10]  = p_x*p_1*cos(p_1*X[n])+p_xy*p_3*Y[n]*cos(p_3*X[n]*Y[n]);
			RUVP[11]  = p_y*p_2*cos(p_2*Y[n])+p_xy*p_3*X[n]*cos(p_3*X[n]*Y[n]);
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}*/

		//MMS2
		/*for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n])+rho_xy*cos(rho_3*X[n]*Y[n]);
			RUVP[1]  =  rho_x*rho_1*cos(rho_1*X[n])-rho_xy*rho_3*Y[n]*sin(rho_3*X[n]*Y[n]);
			RUVP[2]  =  -rho_y*rho_2*sin(rho_2*Y[n])-rho_xy*rho_3*X[n]*sin(rho_3*X[n]*Y[n]);
			RUVP[3]  =  u_0 + u_x*sin(u_1*X[n])+u_y*cos(u_2*Y[n])+u_xy*cos(u_3*X[n]*Y[n]);
			RUVP[4]  =  u_x*u_1*cos(u_1*X[n])-u_xy*u_3*Y[n]*sin(u_3*X[n]*Y[n]);
			RUVP[5]  =  -u_y*u_2*sin(u_2*Y[n])-u_xy*u_3*X[n]*sin(u_3*X[n]*Y[n]);
			RUVP[6]  =  v_0 + v_x*cos(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*cos(v_3*X[n]*Y[n]);
			RUVP[7]  =  -v_x*v_1*sin(v_1*X[n])-v_xy*v_3*Y[n]*sin(v_3*X[n]*Y[n]);
			RUVP[8]  =  v_y*v_2*cos(v_2*Y[n])-v_xy*v_3*X[n]*sin(v_3*X[n]*Y[n]);
			RUVP[9]  =  p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			RUVP[10]  = -p_x*p_1*sin(p_1*X[n])+p_xy*p_3*Y[n]*cos(p_3*X[n]*Y[n]);
			RUVP[11]  = p_y*p_2*cos(p_2*Y[n])+p_xy*p_3*X[n]*cos(p_3*X[n]*Y[n]);
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}*/

		//MMS3(SLIPWALL)
		/*for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			RUVP[1]  =  rho_x*rho_1*cos(rho_1*X[n])+rho_xy*rho_3*Y[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[2]  =  rho_y*rho_2*cos(rho_2*Y[n])+rho_xy*rho_3*X[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[9]  =  p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			RUVP[10]  = p_x*p_1*cos(p_1*X[n])+p_xy*p_3*Y[n]*cos(p_3*X[n]*Y[n]);
			RUVP[11]  = p_y*p_2*cos(p_2*Y[n])+p_xy*p_3*X[n]*cos(p_3*X[n]*Y[n]);

			//Domain transformations and partial derivatives of the transformations.
			double const xi = (2/c)*X[n]-1,
						 eta = (2*a*Y[n])/(b*sqrt(a*a+X[n]*X[n]))-1,
						 xi_x = 2/c,
				   		 eta_x = (-2*a*X[n]*Y[n]/b)*pow(a*a+X[n]*X[n],-1.5),
						 xi_y = 0,
						 eta_y = (2*a)/(b*sqrt(a*a+X[n]*X[n]));

			double const u_a = a*sqrt(a*a+0.25*c*c*(xi+1)*(xi+1)),
						 u_b = u_1*sin(u_2*eta)+u_3*cos(u_4*eta),
						 u_c = f_0 + f_1*sin(f_2*(xi+1))+f_3*sin(f_4*(eta+1))+f_5*sin(f_6*(xi+1)*(eta+1)),
						 v_a = 0.5*b*c*(xi+1),
						 v_b = v_1*sin(v_2*eta)+v_3*cos(v_4*eta),
						 v_c = f_0 + f_1*sin(f_2*(xi+1))+f_3*sin(f_4*(eta+1))+f_5*sin(f_6*(xi+1)*(eta+1)),
						 u_a_xi = 0.25*a*c*c*(xi+1)/sqrt(a*a+0.25*c*c*(xi+1)*(xi+1)),
						 u_a_eta = 0,
					   	 v_a_xi = 0.5*b*c,
						 v_a_eta = 0,
						 u_b_xi = 0,
						 u_b_eta = u_1*u_2*cos(u_2*eta)-u_3*u_4*sin(u_4*eta),
					   	 v_b_xi = 0,
						 v_b_eta = v_1*v_2*cos(v_2*eta)-v_3*v_4*sin(v_4*eta),
						 u_c_xi = f_1*f_2*cos(f_2*(xi+1))+f_5*f_6*(eta+1)*cos(f_6*(xi+1)*(eta+1)),
						 u_c_eta = f_3*f_4*cos(f_4*(eta+1))+f_5*f_6*(xi+1)*cos(f_6*(xi+1)*(eta+1)),
						 v_c_xi = f_1*f_2*cos(f_2*(xi+1))+f_5*f_6*(eta+1)*cos(f_6*(xi+1)*(eta+1)),
						 v_c_eta = f_3*f_4*cos(f_4*(eta+1))+f_5*f_6*(xi+1)*cos(f_6*(xi+1)*(eta+1));

			double const u_a_x = u_a_xi*xi_x+u_a_eta*eta_x,
						 u_a_y = u_a_xi*xi_y+u_a_eta*eta_y,
		 				 u_b_x = u_b_xi*xi_x+u_b_eta*eta_x,
						 u_b_y = u_b_xi*xi_y+u_b_eta*eta_y,
						 u_c_x = u_c_xi*xi_x+u_c_eta*eta_x,
						 u_c_y = u_c_xi*xi_y+u_c_eta*eta_y,
						 v_a_x = v_a_xi*xi_x+v_a_eta*eta_x,
						 v_a_y = v_a_xi*xi_y+v_a_eta*eta_y,
		 				 v_b_x = v_b_xi*xi_x+v_b_eta*eta_x,
						 v_b_y = v_b_xi*xi_y+v_b_eta*eta_y,
						 v_c_x = v_c_xi*xi_x+v_c_eta*eta_x,
						 v_c_y = v_c_xi*xi_y+v_c_eta*eta_y;

			RUVP[3] = u_a*u_b*u_c;
			RUVP[4] = u_a_x*u_b*u_c+u_a*u_b_x*u_c+u_a*u_b*u_c_x;
			RUVP[5] = u_a_y*u_b*u_c+u_a*u_b_y*u_c+u_a*u_b*u_c_y;
			RUVP[6] = v_a*v_b*v_c;
			RUVP[7] = v_a_x*v_b*v_c+v_a*v_b_x*v_c+v_a*v_b*v_c_x;
			RUVP[8] = v_a_y*v_b*v_c+v_a*v_b_y*v_c+v_a*v_b*v_c_y;
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}*/

		//MMS4
		for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_1*cos(rho_2*X[n])*cos(rho_3*Y[n]);
			RUVP[1]  =  -rho_1*rho_2*cos(rho_3*Y[n])*sin(rho_2*X[n]);
			RUVP[2]  =  -rho_1*rho_3*cos(rho_2*X[n])*sin(rho_3*Y[n]);
			RUVP[3]  =  u_0 + u_1*cos(u_2*X[n])*cos(u_3*Y[n]);
			RUVP[4]  =  -u_1*u_2*cos(u_3*Y[n])*sin(u_2*X[n]);
			RUVP[5]  =  -u_1*u_3*cos(u_2*X[n])*sin(u_3*Y[n]);
			RUVP[6]  =  v_0 + v_1*cos(v_2*X[n])*cos(v_3*Y[n]);
			RUVP[7]  =  -v_1*v_2*cos(v_3*Y[n])*sin(v_2*X[n]);
			RUVP[8]  =  -v_1*v_3*cos(v_2*X[n])*sin(v_3*Y[n]);
			RUVP[9]  =  p_0 + p_1*cos(p_2*X[n])*cos(p_3*Y[n]);
			RUVP[10]  =  -p_1*p_2*cos(p_3*Y[n])*sin(p_2*X[n]);
			RUVP[11]  =  -p_1*p_3*cos(p_2*X[n])*sin(p_3*Y[n]);
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}

	} else if (strstr(TestCase,"ParabolicPipe")) {
		if (d != 2)
			EXIT_UNSUPPORTED;

		double rho_0 = DB.rho_store[0],
		       rho_x = DB.rho_store[1],
		       rho_1 = DB.rho_store[2],
		       rho_y = DB.rho_store[3],
		       rho_2 = DB.rho_store[4],
		       rho_xy = DB.rho_store[5],
		       rho_3 = DB.rho_store[6];

		double p_0 = DB.p_store[0],
		       p_x = DB.p_store[1],
		       p_1 = DB.p_store[2],
		       p_y = DB.p_store[3],
		       p_2 = DB.p_store[4],
		       p_xy = DB.p_store[5],
		       p_3 = DB.p_store[6];

		double u_0 = DB.u_store[0],
		       u_x = DB.u_store[1],
		       u_1 = DB.u_store[2],
		       u_y = DB.u_store[3],
		       u_2 = DB.u_store[4],
		       u_xy = DB.u_store[5],
		       u_3 = DB.u_store[6];

	 	double v_0 = DB.v_store[0],
		       v_x = DB.v_store[1],
		       v_1 = DB.v_store[2],
		       v_y = DB.v_store[3],
		       v_2 = DB.v_store[4],
      	       v_xy = DB.v_store[5],
		       v_3 = DB.v_store[6];

		/*double a_1 = DB.geo_store[0],
			   a_2 = DB.geo_store[1],
			   b = DB.geo_store[2];*/

		double RUVP[12];

		/*for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n])+rho_xy*cos(rho_3*X[n]*Y[n]);
			RUVP[1]  =  rho_x*rho_1*cos(rho_1*X[n])-rho_xy*rho_3*Y[n]*sin(rho_3*X[n]*Y[n]);
			RUVP[2]  = -rho_y*rho_2*sin(rho_2*Y[n])-rho_xy*rho_3*X[n]*sin(rho_3*X[n]*Y[n]);
			RUVP[3]  =  u_0 + u_x*sin(u_1*X[n])+u_y*cos(u_2*Y[n])+u_xy*cos(u_3*X[n]*Y[n]);
			RUVP[4]  =  u_x*u_1*cos(u_1*X[n])-u_xy*u_3*Y[n]*sin(u_3*X[n]*Y[n]);
			RUVP[5]  = -u_y*u_2*sin(u_2*Y[n])-u_xy*u_3*X[n]*sin(u_3*X[n]*Y[n]);
			RUVP[6]  =  v_0 + v_x*cos(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*cos(v_3*X[n]*Y[n]);
			RUVP[7]  = -v_x*v_1*sin(v_1*X[n])-v_xy*v_3*Y[n]*sin(v_3*X[n]*Y[n]);
			RUVP[8]  = v_y*v_2*cos(v_2*Y[n])-v_xy*v_3*X[n]*sin(v_3*X[n]*Y[n]);
			RUVP[9]  =  p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			RUVP[10]  = -p_x*p_1*sin(p_1*X[n])+p_xy*p_3*Y[n]*cos(p_3*X[n]*Y[n]);
			RUVP[11]  = p_y*p_2*cos(p_2*Y[n])+p_xy*p_3*X[n]*cos(p_3*X[n]*Y[n]);
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}*/

		for (size_t n = 0; n < Nn; n++) {
			RUVP[0]  =  rho_0 + rho_x*sin(rho_1*X[n])+rho_y*sin(rho_2*Y[n])+rho_xy*sin(rho_3*X[n]*Y[n]);
			RUVP[1]  =  rho_x*rho_1*cos(rho_1*X[n])+rho_xy*rho_3*Y[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[2]  =  rho_y*rho_2*cos(rho_2*Y[n])+rho_xy*rho_3*X[n]*cos(rho_3*X[n]*Y[n]);
			RUVP[3]  =  u_0 + u_x*sin(u_1*X[n])+u_y*sin(u_2*Y[n])+u_xy*sin(u_3*X[n]*Y[n]);
			RUVP[4]  =  u_x*u_1*cos(u_1*X[n])+u_xy*u_3*Y[n]*cos(u_3*X[n]*Y[n]);
			RUVP[5]  =  u_y*u_2*cos(u_2*Y[n])+u_xy*u_3*X[n]*cos(u_3*X[n]*Y[n]);
			RUVP[6]  =  v_0 + v_x*sin(v_1*X[n])+v_y*sin(v_2*Y[n])+v_xy*sin(v_3*X[n]*Y[n]);
			RUVP[7]  =  v_x*v_1*cos(v_1*X[n])+v_xy*v_3*Y[n]*cos(v_3*X[n]*Y[n]);
			RUVP[8]  =  v_y*v_2*cos(v_2*Y[n])+v_xy*v_3*X[n]*cos(v_3*X[n]*Y[n]);
			RUVP[9]  =  p_0 + p_x*sin(p_1*X[n])+p_y*sin(p_2*Y[n])+p_xy*sin(p_3*X[n]*Y[n]);
			RUVP[10]  = p_x*p_1*cos(p_1*X[n])+p_xy*p_3*Y[n]*cos(p_3*X[n]*Y[n]);
			RUVP[11]  = p_y*p_2*cos(p_2*Y[n])+p_xy*p_3*X[n]*cos(p_3*X[n]*Y[n]);
				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}

		/*for (size_t n = 0; n < Nn; n++) {

			double const xi = 2*sqrt(b/a_1)*X[n]-1,
						 eta = (2*Y[n]+2*b*X[n]*X[n]-a_1-a_2)/(a_2-a_1),
						 xi_x = 2*sqrt(b/a_1),
				   		 eta_x = (4*b*X[n])/(a_2-a_1),
						 xi_y = 0,
						 eta_y = 2/(a_2-a_1);

			RUVP[0] = rho_0 + rho_1*sin(rho_2*(xi+rho_3))+rho_4*sin(rho_5*(eta+rho_6));
			RUVP[3] = u_0 + u_1*sin(u_2*(xi+u_3))+u_4*sin(u_5*(eta+u_6));
			RUVP[6] = v_0 + v_1*sin(v_2*(xi+v_3))+v_4*sin(v_5*(eta+v_6));
			RUVP[9] = p_0 + p_1*sin(p_2*(xi+p_3))+p_4*sin(p_5*(eta+p_6));

			double const rho_xi = rho_1*rho_2*cos(rho_2*(xi+rho_3)),
						 rho_eta = rho_4*rho_5*cos(rho_5*(eta+rho_6)),
					   	 u_xi = u_1*u_2*cos(u_2*(xi+u_3)),
						 u_eta = u_4*u_5*cos(u_5*(eta+u_6)),
					   	 v_xi = v_1*v_2*cos(v_2*(xi+v_3)),
						 v_eta = v_4*v_5*cos(v_5*(eta+v_6)),
					   	 p_xi = p_1*p_2*cos(p_2*(xi+p_3)),
						 p_eta = p_4*p_5*cos(p_5*(eta+p_6));

			RUVP[1] = rho_xi*xi_x + rho_eta*eta_x;
			RUVP[2] = rho_xi*xi_y + rho_eta*eta_y;
			RUVP[4] = u_xi*xi_x + u_eta*eta_x;
			RUVP[5] = u_xi*xi_y + u_eta*eta_y;
			RUVP[7] = v_xi*xi_x + v_eta*eta_x;
			RUVP[8] = v_xi*xi_y + v_eta*eta_y;
			RUVP[10] = p_xi*xi_x + p_eta*eta_x;
			RUVP[11] = p_xi*xi_y + p_eta*eta_y;

				for (size_t eq = 0; eq < Neq; eq++)
						source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
		}*/


	} else if (strstr(TestCase,"SinusoidalPipe")) {

             if (d == 2) {
                                double rho_0 = DB.rho_store[0],
                                       rho_x = DB.rho_store[1],
                                       rho_1 = DB.rho_store[2],
                                       rho_y = DB.rho_store[3],
                                       rho_2 = DB.rho_store[4];

                                double p_0 = DB.p_store[0],
                                       p_x = DB.p_store[1],
                                       p_1 = DB.p_store[2],
                                       p_y = DB.p_store[3],
                                       p_2 = DB.p_store[4];

                                double w_0 = DB.w_store[0],
                                       w_x = DB.w_store[1],
                                       w_1 = DB.w_store[2],
                                       w_y = DB.w_store[3],
                                       w_2 = DB.w_store[4];

                                double a = 2,  b = PI/2;

                                double RUVP[12];


                   for (eq = 0; eq < Neq; eq++) {
                           for(n = 0; n < Nn; n++) {
                                    *(RUVP) = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n]);
                                    *(RUVP+1) =  rho_1*rho_x*cos(rho_1*X[n]);
                                    *(RUVP+2) = -rho_2*rho_y*sin(rho_2*Y[n]);
                                    *(RUVP+3) = w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]);
                                    *(RUVP+4) = w_1*w_x*cos(w_1*X[n]);
                                    *(RUVP+5) = -w_2*w_y*sin(w_2*Y[n]);
                                    *(RUVP+6) = (w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]))*(-a*b*sin(b*X[n]));
                                    *(RUVP+7) = -a*b*sin(b*X[n])*(w_1*w_x*cos(w_1*X[n]))-a*pow(b,2)*cos(b*X[n])*(w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]));
                                    *(RUVP+8) = a*b*sin(b*X[n])*w_y*w_2*sin(w_2*Y[n]);
                                    *(RUVP+9) = p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n]);
                                    *(RUVP+10) = -p_1*p_x*sin(p_1*X[n]);
                                    *(RUVP+11) = p_2*p_y*cos(p_2*Y[n]);
                                     source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
                           }
                   }
             } else
                   EXIT_UNSUPPORTED;


       } else {
		EXIT_UNSUPPORTED;
	}
}

static double get_boundary_value_Advection(double const x, double const y, double const z)
{
	char const *const PDESpecifier = DB.PDESpecifier;

	double uEx = 0.0;
	if (strstr(PDESpecifier,"Steady/Default") || strstr(PDESpecifier,"Steady/Peterson")) {
		double const *const b = DB.ADV_b;

		if (!(b[0] == 0.0 && b[1] == 1.0 && b[2] == 0.0))
			EXIT_UNSUPPORTED;

		uEx = sin(2.0*x);
	} else {
		EXIT_UNSUPPORTED;
		printf("%f %f\n",y,z);
	}
	return uEx;
}

static double generate_Euler_source(unsigned int eq, double *RUVP)
{
	if (DB.d != 2)
		EXIT_UNSUPPORTED;

	double r, rx, ry, u, ux, uy, v, vx, vy, p, px, py;
	double f;

	r = RUVP[0], rx = RUVP[1],  ry = RUVP[2];
	u = RUVP[3], ux = RUVP[4],  uy = RUVP[5];
	v = RUVP[6], vx = RUVP[7],  vy = RUVP[8];
	p = RUVP[9], px = RUVP[10], py = RUVP[11];

	if (eq == 1) {
		f = rx*u+ux*r+ry*v+vy*r;
	} else if (eq == 2) {
		f = px+pow(u,2)*rx+2*r*u*ux+r*u*vy+r*v*uy+v*u*ry;
	} else if (eq == 3) {
		f = py+pow(v,2)*ry+2*r*v*vy+r*v*ux+r*u*vx+v*u*rx;
	} else if (eq == 4) {
		double const V2  = u*u+v*v,
		             V2x = 2*(u*ux+v*vx),
		             V2y = 2*(u*uy+v*vy);
		double const E  = p/GM1  + 0.5*r*V2,
		             Ex = px/GM1 + 0.5*(rx*V2+r*V2x),
		             Ey = py/GM1 + 0.5*(ry*V2+r*V2y);
		f = (Ex+px)*u + (E+p)*ux + (Ey+py)*v + (E+p)*vy;
	} else {
		EXIT_UNSUPPORTED;
	}
	return f;
}
