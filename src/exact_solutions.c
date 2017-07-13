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
                           rho_x = DB.rho_store[1],
                           rho_1 = DB.rho_store[2],
                           rho_y = DB.rho_store[3],
                           rho_2 = DB.rho_store[4];

                    double p_0 = DB.p_store[0],
                           p_x = DB.p_store[1],
                           p_1 = DB.p_store[2],
                           p_y = DB.p_store[3],
                           p_2 = DB.p_store[4];

                    double a = 2, b = 4;

                    if (d == 2 ) {
                            for (i = 0; i < Nn; i++) {
                         rhoEx[i] = rho_0 + rho_x*cos(rho_1*X[i])+rho_y*sin(rho_2*Y[i]);
                         pEx[i] = p_0 + p_x*cos(p_1*X[i])+p_y*sin(p_2*Y[i]);
                         uEx[i] = pow(a,2)*Y[i];
                         vEx[i] = -pow(b,2)*X[i];
                         wEx[i] = 0.0;
                            }
                    } else
                           EXIT_UNSUPPORTED;

        } else if (strstr(TestCase,"ParabolicPipe")) {
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

                    double b = 2;

                    if (d == 2 ) {
                                for (i = 0; i < Nn; i++) {
                             rhoEx[i] = rho_0 + rho_x*sin(rho_1*X[i])+rho_y*cos(rho_2*Y[i]);
                             pEx[i] = p_0 + p_x*cos(p_1*X[i])+p_y*sin(p_2*Y[i]);
                             uEx[i] = w_0 + w_x*sin(w_1*X[i])+w_y*cos(w_2*Y[i]);
                             vEx[i] = (-2*b*X[i])*(w_0 + w_x*sin(w_1*X[i])+w_y*cos(w_2*Y[i]));
                             wEx[i] = 0.0;

                                }
                    } else
                            EXIT_UNSUPPORTED;

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
	const double *X, *Y, *Z;

	            X = &XYZ[Nn*0];
                    Y = &XYZ[Nn*1];
                    Z = &XYZ[Nn*(d-1)];


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
			}
		}
	} else if (strstr(TestCase,"EllipticPipe")) {

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

                                double a = 2, b = 4;

                                double RUVP[12];


                   for (eq = 0; eq < Neq; eq++) {
                           for(n = 0; n < Nn; n++) {
                                    *(RUVP) = rho_0 + rho_x*cos(rho_1*X[n])+rho_y*sin(rho_2*Y[n]);
                                    *(RUVP+1) =  -rho_1*rho_x*sin(rho_1*X[n]);
                                    *(RUVP+2) = rho_2*rho_y*cos(rho_2*Y[n]);
                                    *(RUVP+3) = pow(a,2)*Y[n];
                                    *(RUVP+4) = 0;
                                    *(RUVP+5) = pow(a,2);
                                    *(RUVP+6) = -pow(b,2)*X[n];
                                    *(RUVP+7) = -pow(b,2);
                                    *(RUVP+8) = 0;
                                    *(RUVP+9) = p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n]);
                                    *(RUVP+10) = -p_1*p_x*sin(p_1*X[n]);
                                    *(RUVP+11) = p_2*p_y*cos(p_2*Y[n]);
                                     source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
                           }
                   }
            						printf("error\n");
			} else
                   EXIT_UNSUPPORTED;


       } else if (strstr(TestCase,"ParabolicPipe")) {

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

                                double  b = 2;

                                double RUVP[12];


                   for (eq = 0; eq < Neq; eq++) {
                           for(n = 0; n < Nn; n++) {
                                    *(RUVP) = rho_0 + rho_x*sin(rho_1*X[n])+rho_y*cos(rho_2*Y[n]);
                                    *(RUVP+1) =  rho_1*rho_x*cos(rho_1*X[n]);
                                    *(RUVP+2) = -rho_2*rho_y*sin(rho_2*Y[n]);
                                    *(RUVP+3) = w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]);
                                    *(RUVP+4) = w_1*w_x*cos(w_1*X[n]);
                                    *(RUVP+5) = -w_2*w_y*sin(w_2*Y[n]);
                                    *(RUVP+6) = (w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]))*(-2*b*X[n]);
                                    *(RUVP+7) = -2*b*X[n]*(w_1*w_x*cos(w_1*X[n]))-2*b*(w_0+w_x*sin(w_1*X[n])+w_y*cos(w_2*Y[n]));
                                    *(RUVP+8) = 2*b*X[n]*w_y*w_2*sin(w_2*Y[n]);
                                    *(RUVP+9) = p_0 + p_x*cos(p_1*X[n])+p_y*sin(p_2*Y[n]);
                                    *(RUVP+10) = -p_1*p_x*sin(p_1*X[n]);
                                    *(RUVP+11) = p_2*p_y*cos(p_2*Y[n]);
                                     source[eq*Nn+n] = generate_Euler_source(eq+1, RUVP);
                           }
                   }
             } else
                   EXIT_UNSUPPORTED;


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

double r, rx, ry, u, ux, uy, v, vx, vy, p, px, py;
double f;

        r = RUVP[0], rx = RUVP[1], ry = RUVP[2];
        u = RUVP[3], ux = RUVP[4], uy = RUVP[5];
        v = RUVP[6], vx = RUVP[7], vy = RUVP[8];
        p = RUVP[9], px = RUVP[10], py = RUVP[11];

          if (eq == 1) {
           f = rx*u+ux*r+ry*v+vy*r;

        } else if (eq == 2) {
           f = px+pow(u,2)*rx+2*r*u*ux+r*u*vy+r*v*uy+v*u*ry;

        } else if (eq == 3) {
           f = py+pow(v,2)*ry+2*r*v*vy+r*v*ux+r*u*vx+v*u*rx;

        } else if (eq == 4) {
           f = (GAMMA/(GAMMA-1))*(px*u+ux*p)+rx*0.5*(pow(u,3)+u*pow(v,2))+0.5*r*(3*pow(u,2)*ux+pow(v,2)*ux+2*v*u*vx)+
               (GAMMA/(GAMMA-1))*(py*v+vy*p)+ry*0.5*(pow(v,3)+v*pow(u,2))+0.5*r*(3*pow(v,2)*vy+pow(u,2)*vy+2*v*u*uy);
        } else {
           EXIT_UNSUPPORTED; }

        return f;
}
