#include <stdlib.h>
#include <stdio.h>

#include "mkl.h"
/*
 *	Purpose:
 *		Provide functions for testing mm_CTN.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void mm_CTN_mv1_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int mMax, kMax;
	register double *pA = A, *pB = B, *pC = C;

	// First row of A
	*pC = (*pA)*(*pB);
	for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
		pA++;
		pB++;

		*pC += (*pA)*(*pB);
	}

	// Remaining rows of A
	for (mMax = m-1; mMax--; ) { // loop over rows of A
		pA++;
		pC++;
		pB = B;

		*pC = (*pA)*(*pB);
		for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
			pA++;
			pB++;

			*pC += (*pA)*(*pB);
		}
	}
}

void mm_CTN_mm1_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int mMax, nMax, kMax;
	register double *pA = A, *pB = B, *pC = C;

	// First column of B/C
	// First row of A
	*pC = (*pA)*(*pB);
	for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
		pA++;
		pB++;

		*pC += (*pA)*(*pB);
	}

	// Remaining rows of A
	for (mMax = m-1; mMax--; ) { // loop over rows of A
		pA++;
		pC++;
		pB = B;

		*pC = (*pA)*(*pB);
		for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
			pA++;
			pB++;

			*pC += (*pA)*(*pB);
		}
	}

	// Remaining columns of B/C
	for (nMax = n-1; nMax--; ) { // loop over columns of B/C
		pA = A;
		pB = B+(n-nMax-1)*k;
		pC++;
		*pC = (*pA)*(*pB);

		for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
			pA++;
			pB++;

			*pC += (*pA)*(*pB);
		}

		// Remaining rows of A
		for (mMax = m-1; mMax--; ) { // loop over rows of A
			pA++;
			pC++;
			pB = B+(n-nMax-1)*k;

			*pC = (*pA)*(*pB);
			for (kMax = k-1; kMax--; ) { // loop over columns of A/rows of B/C
				pA++;
				pB++;

				*pC += (*pA)*(*pB);
			}
		}
	}
}

void mm_CTN_mv2_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int i, j, iMax, jMax, IndA;
	for (i = 0, iMax = m; i < iMax; i++) {
		C[i] = 0.0;
		IndA = i*k;
		for (j = 0, jMax = k; j < jMax; j++) {
			C[i] += A[IndA+j]*B[j];
		}
	}
}

void mm_CTN_mm2_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int i, j, l, iMax, jMax, lMax, IndA, IndB, IndC;

	for (i = 0, iMax = n; i < iMax; i++) {
	for (j = 0, jMax = m; j < jMax; j++) {
		IndC = i*m+j;
		IndA = j*k;
		IndB = i*k;

		C[IndC] = 0.0;
		for (l = 0, lMax = k; l < lMax; l++) {
			C[IndC] += A[IndA+l]*B[IndB+l];
		}
	}}
}

void mm_CTN_mv3_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int i, j, iMax, jMax, IndA;
	for (i = 0, iMax = m; i < iMax; i++) {
		IndA = i*k;
		C[i] = A[IndA]*B[0];
		for (j = 1, jMax = k; j < jMax; j++) {
			C[i] += A[IndA+j]*B[j];
		}
	}
}

void mm_CTN_mm3_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int i, j, l, iMax, jMax, lMax, IndA, IndB, IndC;

	for (i = 0, iMax = n; i < iMax; i++) {
	for (j = 0, jMax = m; j < jMax; j++) {
		IndC = i*m+j;
		IndA = j*k;
		IndB = i*k;

		C[IndC] = A[IndA]*B[IndB];
		for (l = 1, lMax = k; l < lMax; l++) {
			C[IndC] += A[IndA+l]*B[IndB+l];
		}
	}}
}

void mm_CTN_mvBLAS_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	MKL_INT m_MKL   = k,
			k_MKL   = m,
			inc_MKL = 1;

	cblas_dgemv(CblasColMajor,CblasTrans,m_MKL,k_MKL,1.0,A,m_MKL,B,inc_MKL,0.0,C,inc_MKL);
}

void mm_CTN_mmBLAS_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	MKL_INT m_MKL = m,
			n_MKL = n,
			k_MKL = k;

	cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,m_MKL,n_MKL,k_MKL,1.0,A,k_MKL,B,k_MKL,0.0,C,m_MKL);
}

void mm_CTN_mv4_5d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	register unsigned int kMax;
	register double *a0 = A, *a1 = A+5, *a2 = A+10, *a3 = A+15, *a4 = A+20;
	register double *b = B;
	register double *c0 = C, *c1 = C+1, *c2 = C+2, *c3 = C+3, *c4 = C+4;

	// Loop unrolled version of mv multiplication for m = k = 5, n = 1

	// First column
	*c0 = (*a0)*(*b);
	*c1 = (*a1)*(*b);
	*c2 = (*a2)*(*b);
	*c3 = (*a3)*(*b);
	*c4 = (*a4)*(*b);

	// Remaining columns
	for (kMax = k-1; kMax--; ) { // loop over columns of A
		a0++;
		a1++;
		a2++;
		a3++;
		a4++;
		b++;

		*c0 += (*a0)*(*b);
		*c1 += (*a1)*(*b);
		*c2 += (*a2)*(*b);
		*c3 += (*a3)*(*b);
		*c4 += (*a4)*(*b);
	}
}

void mm_CTN_mv_unrolled_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	// Square A inputs only (for the time being)!

	register unsigned int kMax;

	switch(k) {
		case 1: {
			break;
		} case 2: {
			break;
		} case 3: {
			break;
		} case 4: {
			break;
		} case 5: {
			register double *a0 = A, *a1 = A+5, *a2 = A+10, *a3 = A+15, *a4 = A+20;
			register double *b = B;
			register double *c0 = C, *c1 = C+1, *c2 = C+2, *c3 = C+3, *c4 = C+4;

			// First column
			*c0 = (*a0)*(*b), *c1 = (*a1)*(*b), *c2 = (*a2)*(*b), *c3 = (*a3)*(*b), *c4 = (*a4)*(*b);

			// Remaining columns
			for (kMax = k-1; kMax--; ) { // loop over columns of A
				a0++, a1++, a2++, a3++, a4++;
				b++;

				*c0 += (*a0)*(*b), *c1 += (*a1)*(*b), *c2 += (*a2)*(*b), *c3 += (*a3)*(*b), *c4 += (*a4)*(*b);
			}
			break;
		} case 6: {
			break;
		} case 7: {
			break;
		} case 8: {
			register double *a0 = A, *a1 = A+8, *a2 = A+16, *a3 = A+24, *a4 = A+32, *a5 = A+40, *a6 = A+48, *a7 = A+56;
			register double *b = B;
			register double *c0 = C, *c1 = C+1, *c2 = C+2, *c3 = C+3, *c4 = C+4, *c5 = C+5, *c6 = C+6, *c7 = C+7;

			// First column
			*c0 = (*a0)*(*b), *c1 = (*a1)*(*b), *c2 = (*a2)*(*b), *c3 = (*a3)*(*b),
			*c4 = (*a4)*(*b), *c5 = (*a5)*(*b), *c6 = (*a6)*(*b), *c7 = (*a7)*(*b);

			// Remaining columns
			for (kMax = k-1; kMax--; ) { // loop over columns of A
				a0++, a1++, a2++, a3++, a4++, a5++, a6++, a7++;
				b++;

				*c0 += (*a0)*(*b), *c1 += (*a1)*(*b), *c2 += (*a2)*(*b), *c3 += (*a3)*(*b),
				*c4 += (*a4)*(*b), *c5 += (*a5)*(*b), *c6 += (*a6)*(*b), *c7 += (*a7)*(*b);
			}
			break;
		} default: {
			printf("Error: unsupported size in mm_CTN_mv_unrolled_d.\n"), exit(1);
			break;
		}
	}
}

void mm_CTN_mv_fully_unrolled_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	// Square A inputs only (for the time being)!
	// C/A values only used once, b values used many times => for mm, A values used many times (generalize here)

	switch(m) {
	case 1: {
		break;
	} case 2: {
		break;
	} case 3: {
		break;
	} case 4: {
		switch(k) {
		case 1: {
			break;
		} case 4: {
			register double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
							*a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
							*a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
							*a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15;
			register double *b0 = B    , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ;

			C[0] = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3);
			C[1] = (*a4 )*(*b0) + (*a5 )*(*b1) + (*a6 )*(*b2) + (*a7 )*(*b3);
			C[2] = (*a8 )*(*b0) + (*a9 )*(*b1) + (*a10)*(*b2) + (*a11)*(*b3);
			C[3] = (*a12)*(*b0) + (*a13)*(*b1) + (*a14)*(*b2) + (*a15)*(*b3);

			break;
		} default: {
			printf("message m4\n"), exit(1);
			break;
		}}
		break;
	} case 5: {
		switch(k) {
		case 1: {
			break;
		} case 2: {
			break;
		} case 3: {
			break;
		} case 4: {
			break;
		} case 5: {
			register double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
							*a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
							*a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
							*a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
							*a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24;
			register double *b0 = B    , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;

			C[0] = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4);
			C[1] = (*a5 )*(*b0) + (*a6 )*(*b1) + (*a7 )*(*b2) + (*a8 )*(*b3) + (*a9 )*(*b4);
			C[2] = (*a10)*(*b0) + (*a11)*(*b1) + (*a12)*(*b2) + (*a13)*(*b3) + (*a14)*(*b4);
			C[3] = (*a15)*(*b0) + (*a16)*(*b1) + (*a17)*(*b2) + (*a18)*(*b3) + (*a19)*(*b4);
			C[4] = (*a20)*(*b0) + (*a21)*(*b1) + (*a22)*(*b2) + (*a23)*(*b3) + (*a24)*(*b4);

			break;
		} default: {
			printf("message m5\n"), exit(1);
			break;
		}}
		break;
	} case 8: {
		switch(k) {
		case 1: {
			break;
		} case 2: {
			break;
		} case 3: {
			break;
		} case 4: {
			break;
		} case 5: {
			break;
		} case 6: {
			break;
		} case 7: {
			break;
		} case 8: {
			register double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
							*a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
							*a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
							*a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
							*a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
							*a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
							*a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
							*a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
							*a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35,
							*a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
							*a40 = A+40, *a41 = A+41, *a42 = A+42, *a43 = A+43,
							*a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
							*a48 = A+48, *a49 = A+49, *a50 = A+50, *a51 = A+51,
							*a52 = A+52, *a53 = A+53, *a54 = A+54, *a55 = A+55,
							*a56 = A+56, *a57 = A+57, *a58 = A+58, *a59 = A+59,
							*a60 = A+60, *a61 = A+61, *a62 = A+62, *a63 = A+63;
			register double *b0 = B    , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ,
							*b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;

			C[0] = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4) + (*a5 )*(*b5) + (*a6 )*(*b6) + (*a7 )*(*b7);
			C[1] = (*a8 )*(*b0) + (*a9 )*(*b1) + (*a10)*(*b2) + (*a11)*(*b3) + (*a12)*(*b4) + (*a13)*(*b5) + (*a14)*(*b6) + (*a15)*(*b7);
			C[2] = (*a16)*(*b0) + (*a17)*(*b1) + (*a18)*(*b2) + (*a19)*(*b3) + (*a20)*(*b4) + (*a21)*(*b5) + (*a22)*(*b6) + (*a23)*(*b7);
			C[3] = (*a24)*(*b0) + (*a25)*(*b1) + (*a26)*(*b2) + (*a27)*(*b3) + (*a28)*(*b4) + (*a29)*(*b5) + (*a30)*(*b6) + (*a31)*(*b7);
			C[4] = (*a32)*(*b0) + (*a33)*(*b1) + (*a34)*(*b2) + (*a35)*(*b3) + (*a36)*(*b4) + (*a37)*(*b5) + (*a38)*(*b6) + (*a39)*(*b7);
			C[5] = (*a40)*(*b0) + (*a41)*(*b1) + (*a42)*(*b2) + (*a43)*(*b3) + (*a44)*(*b4) + (*a45)*(*b5) + (*a46)*(*b6) + (*a47)*(*b7);
			C[6] = (*a48)*(*b0) + (*a49)*(*b1) + (*a50)*(*b2) + (*a51)*(*b3) + (*a52)*(*b4) + (*a53)*(*b5) + (*a54)*(*b6) + (*a55)*(*b7);
			C[7] = (*a56)*(*b0) + (*a57)*(*b1) + (*a58)*(*b2) + (*a59)*(*b3) + (*a60)*(*b4) + (*a61)*(*b5) + (*a62)*(*b6) + (*a63)*(*b7);

			break;
		} default: {
			printf("Error: unsupported size in mm_CTN_mv_unrolled_d.\n"), exit(1);
			break;
		}}
		break;
	} default: {
		printf("Error: unsupported size in mm_CTN_mv_unrolled_d.\n"), exit(1);
		break;
	}}
}

void mm_CTN_fully_unrolled_mv_d(const int m, const int n, const int k, double *A, double *B, double *C)
{
	/*
	 *	Fully unrolled matrix-vector operations for supported cases.
	 *	Nearly identical timing to mv_fully_unrolled.
	 *	When profiling, figure out which cases are used most frequently and compare with BLAS (ToBeDeleted).
	 */

	switch(m) {
	case 1: {
		break; // m1
	} case 2: {
		break; // m2
	} case 3: {
		break; // m3
	} case 4: {
		break; // m4
	} case 5: {
		switch(k) {
		case 1: {
			break; // m5k1
		} case 2: {
			break; // m5k2
		} case 3: {
			break; // m5k3
		} case 4: {
			break; // m5k4
		} case 5: {
			register double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 , *a4  = A+4 ,
							*a5  = A+5 , *a6  = A+6 , *a7  = A+7 , *a8  = A+8 , *a9  = A+9 ,
							*a10 = A+10, *a11 = A+11, *a12 = A+12, *a13 = A+13, *a14 = A+14,
							*a15 = A+15, *a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
							*a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23, *a24 = A+24;
			register double *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 , *b4  = B+4 ;
			register unsigned int nRem = n-1;

			*C     = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4);
			*(++C) = (*a5 )*(*b0) + (*a6 )*(*b1) + (*a7 )*(*b2) + (*a8 )*(*b3) + (*a9 )*(*b4);
			*(++C) = (*a10)*(*b0) + (*a11)*(*b1) + (*a12)*(*b2) + (*a13)*(*b3) + (*a14)*(*b4);
			*(++C) = (*a15)*(*b0) + (*a16)*(*b1) + (*a17)*(*b2) + (*a18)*(*b3) + (*a19)*(*b4);
			*(++C) = (*a20)*(*b0) + (*a21)*(*b1) + (*a22)*(*b2) + (*a23)*(*b3) + (*a24)*(*b4);

			for ( ; nRem--; ) {
				b0 += k, b1 += k, b2 += k, b3 += k, b4 += k;

				*(++C) = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4);
				*(++C) = (*a5 )*(*b0) + (*a6 )*(*b1) + (*a7 )*(*b2) + (*a8 )*(*b3) + (*a9 )*(*b4);
				*(++C) = (*a10)*(*b0) + (*a11)*(*b1) + (*a12)*(*b2) + (*a13)*(*b3) + (*a14)*(*b4);
				*(++C) = (*a15)*(*b0) + (*a16)*(*b1) + (*a17)*(*b2) + (*a18)*(*b3) + (*a19)*(*b4);
				*(++C) = (*a20)*(*b0) + (*a21)*(*b1) + (*a22)*(*b2) + (*a23)*(*b3) + (*a24)*(*b4);
			}
			break; // m5k5
		} default: {
			printf("Error: Unsupported m = %d, k = %d in mm_CTN_fully_unrolled.\n",m,k), exit(1);
			break; // m5kd
		}}
		break; // m5
	} case 6: {
		break; // m6
	} case 7: {
		break; // m7
	} case 8: {
		switch(k) {
		case 1: {
			break; // m8k1
		} case 2: {
			break; // m8k2
		} case 3: {
			break; // m8k3
		} case 4: {
			break; // m8k4
		} case 5: {
			break; // m8k5
		} case 6: {
			break; // m8k6
		} case 7: {
			break; // m8k7
		} case 8: {
			register double *a0  = A   , *a1  = A+1 , *a2  = A+2 , *a3  = A+3 ,
							*a4  = A+4 , *a5  = A+5 , *a6  = A+6 , *a7  = A+7 ,
							*a8  = A+8 , *a9  = A+9 , *a10 = A+10, *a11 = A+11,
							*a12 = A+12, *a13 = A+13, *a14 = A+14, *a15 = A+15,
							*a16 = A+16, *a17 = A+17, *a18 = A+18, *a19 = A+19,
							*a20 = A+20, *a21 = A+21, *a22 = A+22, *a23 = A+23,
							*a24 = A+24, *a25 = A+25, *a26 = A+26, *a27 = A+27,
							*a28 = A+28, *a29 = A+29, *a30 = A+30, *a31 = A+31,
							*a32 = A+32, *a33 = A+33, *a34 = A+34, *a35 = A+35,
							*a36 = A+36, *a37 = A+37, *a38 = A+38, *a39 = A+39,
							*a40 = A+40, *a41 = A+41, *a42 = A+42, *a43 = A+43,
							*a44 = A+44, *a45 = A+45, *a46 = A+46, *a47 = A+47,
							*a48 = A+48, *a49 = A+49, *a50 = A+50, *a51 = A+51,
							*a52 = A+52, *a53 = A+53, *a54 = A+54, *a55 = A+55,
							*a56 = A+56, *a57 = A+57, *a58 = A+58, *a59 = A+59,
							*a60 = A+60, *a61 = A+61, *a62 = A+62, *a63 = A+63;
			register double *b0  = B   , *b1  = B+1 , *b2  = B+2 , *b3  = B+3 ,
							*b4  = B+4 , *b5  = B+5 , *b6  = B+6 , *b7  = B+7 ;
			register unsigned int nRem = n-1;

			*C     = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4) + (*a5 )*(*b5) + (*a6 )*(*b6) + (*a7 )*(*b7);
			*(++C) = (*a8 )*(*b0) + (*a9 )*(*b1) + (*a10)*(*b2) + (*a11)*(*b3) + (*a12)*(*b4) + (*a13)*(*b5) + (*a14)*(*b6) + (*a15)*(*b7);
			*(++C) = (*a16)*(*b0) + (*a17)*(*b1) + (*a18)*(*b2) + (*a19)*(*b3) + (*a20)*(*b4) + (*a21)*(*b5) + (*a22)*(*b6) + (*a23)*(*b7);
			*(++C) = (*a24)*(*b0) + (*a25)*(*b1) + (*a26)*(*b2) + (*a27)*(*b3) + (*a28)*(*b4) + (*a29)*(*b5) + (*a30)*(*b6) + (*a31)*(*b7);
			*(++C) = (*a32)*(*b0) + (*a33)*(*b1) + (*a34)*(*b2) + (*a35)*(*b3) + (*a36)*(*b4) + (*a37)*(*b5) + (*a38)*(*b6) + (*a39)*(*b7);
			*(++C) = (*a40)*(*b0) + (*a41)*(*b1) + (*a42)*(*b2) + (*a43)*(*b3) + (*a44)*(*b4) + (*a45)*(*b5) + (*a46)*(*b6) + (*a47)*(*b7);
			*(++C) = (*a48)*(*b0) + (*a49)*(*b1) + (*a50)*(*b2) + (*a51)*(*b3) + (*a52)*(*b4) + (*a53)*(*b5) + (*a54)*(*b6) + (*a55)*(*b7);
			*(++C) = (*a56)*(*b0) + (*a57)*(*b1) + (*a58)*(*b2) + (*a59)*(*b3) + (*a60)*(*b4) + (*a61)*(*b5) + (*a62)*(*b6) + (*a63)*(*b7);

			for ( ; nRem--; ) {
				b0 += k, b1 += k, b2 += k, b3 += k, b4 += k, b5 += k, b6 += k, b7 += k;

				*(++C) = (*a0 )*(*b0) + (*a1 )*(*b1) + (*a2 )*(*b2) + (*a3 )*(*b3) + (*a4 )*(*b4) + (*a5 )*(*b5) + (*a6 )*(*b6) + (*a7 )*(*b7);
				*(++C) = (*a8 )*(*b0) + (*a9 )*(*b1) + (*a10)*(*b2) + (*a11)*(*b3) + (*a12)*(*b4) + (*a13)*(*b5) + (*a14)*(*b6) + (*a15)*(*b7);
				*(++C) = (*a16)*(*b0) + (*a17)*(*b1) + (*a18)*(*b2) + (*a19)*(*b3) + (*a20)*(*b4) + (*a21)*(*b5) + (*a22)*(*b6) + (*a23)*(*b7);
				*(++C) = (*a24)*(*b0) + (*a25)*(*b1) + (*a26)*(*b2) + (*a27)*(*b3) + (*a28)*(*b4) + (*a29)*(*b5) + (*a30)*(*b6) + (*a31)*(*b7);
				*(++C) = (*a32)*(*b0) + (*a33)*(*b1) + (*a34)*(*b2) + (*a35)*(*b3) + (*a36)*(*b4) + (*a37)*(*b5) + (*a38)*(*b6) + (*a39)*(*b7);
				*(++C) = (*a40)*(*b0) + (*a41)*(*b1) + (*a42)*(*b2) + (*a43)*(*b3) + (*a44)*(*b4) + (*a45)*(*b5) + (*a46)*(*b6) + (*a47)*(*b7);
				*(++C) = (*a48)*(*b0) + (*a49)*(*b1) + (*a50)*(*b2) + (*a51)*(*b3) + (*a52)*(*b4) + (*a53)*(*b5) + (*a54)*(*b6) + (*a55)*(*b7);
				*(++C) = (*a56)*(*b0) + (*a57)*(*b1) + (*a58)*(*b2) + (*a59)*(*b3) + (*a60)*(*b4) + (*a61)*(*b5) + (*a62)*(*b6) + (*a63)*(*b7);
			}
			break; // m8k8
		} default: {
			printf("Error: Unsupported m = %d, k = %d in mm_CTN_fully_unrolled.\n",m,k), exit(1);
			break; // m8kd
		}}
		break; // m8
	} default: {
		printf("Error: Unsupported m = %d, k = %d in mm_CTN_fully_unrolled.\n",m,k), exit(1);
		break; // md
	}}
}
