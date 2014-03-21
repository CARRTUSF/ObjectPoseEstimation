#ifndef __EIGENSOLVER_H
#define __EIGENSOLVER_H

/*///////////////////////////////////////////////////////////////////////////////////////
//
//                      INTEL CORPORATION PROPRIETARY INFORMATION
//         This software is supplied under the terms of a license agreement or
//         nondisclosure agreement with Intel Corporation and may not be copied
//         or disclosed except in accordance with the terms of that agreement.
//               Copyright (c) 1999 Intel Corporation. All Rights Reserved.
//
//    RCS:
//       $Source$
//       $Revision: 208 $
//      Purpose:
//      Contents:
//      Authors:
//        Dmitry Abrosimov
//*/

#include <math.h>

/*///////////////////////////////////////////////////////////////////////////////////////
//    Names:      cvJacobiEigens_32f, cvJacobiEigens_64d
//    Purpose:    Eigenvalues & eigenvectors calculation of a symmetric matrix:
//                A Vi  =  Ei Vi
//    Context:
//    Parameters: A(n, n) - source symmetric matrix (n - rows & columns number),
//                V(n, n) - matrix of its eigenvectors
//                          (i-th row is an eigenvector Vi),
//                E(n)    - vector of its eigenvalues
//                          (i-th element is an eigenvalue Ei),
//                eps     - accuracy of diagonalization.
//
//    Returns:
//    CV_NO_ERROR or error code
//    Notes:
//        1. The functions destroy source matrix A, so if you need it after
//           running cvJacobiEigens_???, you must copy it before.
//        2. Eigenvalies and eigenvectors are sorted in Ei absolute value descending.
//        3. Calculation time depends on eps value. If the time isn't very important,
//           we recommend to set eps = 0.
//*/

/*=========================== Single precision function ================================*/

int cvJacobiEigens_32f(float* A, float* V, float* E, int n, float eps) {
	int i, j, k, ind;
	float *AA = A, *VV = V;
	double Amax, anorm = 0, ax;

	/*if ( A == NULL || V == NULL || E == NULL ) return CV_NULLPTR_ERR;*/
	/*if ( n <= 0 )                              return CV_BADSIZE_ERR;*/
	if (eps < 1.0e-7f)
		eps = 1.0e-7f;

	/*-------- Prepare --------*/
	for (i = 0; i < n; i++, VV += n, AA += n) {
		for (j = 0; j < i; j++) {
			double Am = AA[j];
			anorm += Am * Am;
		}
		for (j = 0; j < n; j++)
			VV[j] = 0.f;
		VV[i] = 1.f;
	}

	anorm = sqrt(anorm + anorm);
	ax = anorm * eps / n;
	Amax = anorm;

	while (Amax > ax) {
		Amax /= n;
		do /* while (ind) */
		{
			int p, q;
			float *V1 = V, *A1 = A;
			ind = 0;
			for (p = 0; p < n - 1; p++, A1 += n, V1 += n) {
				float * A2 = A + n * (p + 1), *V2 = V + n * (p + 1);
				for (q = p + 1; q < n; q++, A2 += n, V2 += n) {
					double x, y, c, s, c2, s2, a;
					float *A3, Apq = A1[q], App, Aqq, Aip, Aiq, Vpi, Vqi;
					if (fabs(Apq) < Amax)
						continue;

					ind = 1;

					/*---- Calculation of rotation angle's sine & cosine ----*/
					App = A1[p];
					Aqq = A2[q];
					y = 5.0e-1 * (App - Aqq);
					x = -Apq / sqrt(Apq * Apq + y * y);
					if (y < 0.0)
						x = -x;
					s = x / sqrt(2.0 * (1.0 + sqrt(1.0 - x * x)));
					s2 = s * s;
					c = sqrt(1.0 - s2);
					c2 = c * c;
					a = 2.0 * Apq * c * s;

					/*---- Apq annulation ----*/
					A3 = A;
					for (i = 0; i < p; i++, A3 += n) {
						Aip = A3[p];
						Aiq = A3[q];
						Vpi = V1[i];
						Vqi = V2[i];
						A3[p] = (float) (Aip * c - Aiq * s);
						A3[q] = (float) (Aiq * c + Aip * s);
						V1[i] = (float) (Vpi * c - Vqi * s);
						V2[i] = (float) (Vqi * c + Vpi * s);
					}
					for (; i < q; i++, A3 += n) {
						Aip = A1[i];
						Aiq = A3[q];
						Vpi = V1[i];
						Vqi = V2[i];
						A1[i] = (float) (Aip * c - Aiq * s);
						A3[q] = (float) (Aiq * c + Aip * s);
						V1[i] = (float) (Vpi * c - Vqi * s);
						V2[i] = (float) (Vqi * c + Vpi * s);
					}
					for (; i < n; i++) {
						Aip = A1[i];
						Aiq = A2[i];
						Vpi = V1[i];
						Vqi = V2[i];
						A1[i] = (float) (Aip * c - Aiq * s);
						A2[i] = (float) (Aiq * c + Aip * s);
						V1[i] = (float) (Vpi * c - Vqi * s);
						V2[i] = (float) (Vqi * c + Vpi * s);
					}
					A1[p] = (float) (App * c2 + Aqq * s2 - a);
					A2[q] = (float) (App * s2 + Aqq * c2 + a);
					A1[q] = A2[p] = 0.0f;
				} /*q*/
			} /*p*/
		} while (ind);
		Amax /= n;
	} /* while ( Amax > ax ) */

	for (i = 0, k = 0; i < n; i++, k += n + 1)
		E[i] = A[k];
	/*printf(" M = %d\n", M);*/

	/* -------- ordering --------*/
	for (i = 0; i < n; i++) {
		int m = i;
		float Em = (float) fabs(E[i]);
		for (j = i + 1; j < n; j++) {
			float Ej = (float) fabs(E[j]);
			m = (Em < Ej) ? j : m;
			Em = (Em < Ej) ? Ej : Em;
		}
		if (m != i) {
			int l;
			float b = E[i];
			E[i] = E[m];
			E[m] = b;
			for (j = 0, k = i * n, l = m * n; j < n; j++, k++, l++) {
				b = V[k];
				V[k] = V[l];
				V[l] = b;
			}
		}
	}

	return 0;
}

/*=========================== Double precision function ================================*/

int cvJacobiEigens_64d(double* A, double* V, double* E, int n, double eps) {
	int i, j, k, p, q, ind;
	double *A1 = A, *V1 = V, *A2 = A, *V2 = V;
	double Amax = 0.0, anorm = 0.0, ax, deps;

	/*if ( A == NULL || V == NULL || E == NULL ) return CV_NULLPTR_ERR;*/
	/*if ( n <= 0 )                              return CV_BADSIZE_ERR;*/
	if (eps < 1.0e-15)
		eps = 1.0e-15;
	deps = eps / (double) n;

	/*-------- Prepare --------*/
	for (i = 0; i < n; i++, V1 += n, A1 += n) {
		for (j = 0; j < i; j++) {
			double Am = A1[j];
			anorm += Am * Am;
		}
		for (j = 0; j < n; j++)
			V1[j] = 0.0;
		V1[i] = 1.0;
	}

	anorm = sqrt(anorm + anorm);
	ax = anorm * eps / n;
	Amax = anorm;

	while (Amax > ax) {
		Amax /= n;
		do /* while (ind) */
		{
			ind = 0;
			A1 = A;
			V1 = V;
			for (p = 0; p < n - 1; p++, A1 += n, V1 += n) {
				A2 = A + n * (p + 1);
				V2 = V + n * (p + 1);
				for (q = p + 1; q < n; q++, A2 += n, V2 += n) {
					double x, y, c, s, c2, s2, a;
					double *A3, Apq, App, Aqq, App2, Aqq2, Aip, Aiq, Vpi, Vqi;
					if (fabs(A1[q]) < Amax)
						continue;
					Apq = A1[q];

					ind = 1;

					/*---- Calculation of rotation angle's sine & cosine ----*/
					App = A1[p];
					Aqq = A2[q];
					y = 5.0e-1 * (App - Aqq);
					x = -Apq / sqrt(Apq * Apq + y * y);
					if (y < 0.0)
						x = -x;
					s = x / sqrt(2.0 * (1.0 + sqrt(1.0 - x * x)));
					s2 = s * s;
					c = sqrt(1.0 - s2);
					c2 = c * c;
					a = 2.0 * Apq * c * s;

					/*---- Apq annulation ----*/
					A3 = A;
					for (i = 0; i < p; i++, A3 += n) {
						Aip = A3[p];
						Aiq = A3[q];
						Vpi = V1[i];
						Vqi = V2[i];
						A3[p] = Aip * c - Aiq * s;
						A3[q] = Aiq * c + Aip * s;
						V1[i] = Vpi * c - Vqi * s;
						V2[i] = Vqi * c + Vpi * s;
					}
					for (; i < q; i++, A3 += n) {
						Aip = A1[i];
						Aiq = A3[q];
						Vpi = V1[i];
						Vqi = V2[i];
						A1[i] = Aip * c - Aiq * s;
						A3[q] = Aiq * c + Aip * s;
						V1[i] = Vpi * c - Vqi * s;
						V2[i] = Vqi * c + Vpi * s;
					}
					for (; i < n; i++) {
						Aip = A1[i];
						Aiq = A2[i];
						Vpi = V1[i];
						Vqi = V2[i];
						A1[i] = Aip * c - Aiq * s;
						A2[i] = Aiq * c + Aip * s;
						V1[i] = Vpi * c - Vqi * s;
						V2[i] = Vqi * c + Vpi * s;
					}
					App2 = App * c2 + Aqq * s2 - a;
					Aqq2 = App * s2 + Aqq * c2 + a;
					A1[p] = App2;
					A2[q] = Aqq2;
					A1[q] = A2[p] = 0.0;
				} /*q*/
			} /*p*/
		} while (ind);
	} /* while ( Amax > ax ) */

	for (i = 0, k = 0; i < n; i++, k += n + 1)
		E[i] = A[k];

	/* -------- ordering --------*/
	/*
	for (i = 0; i < n; i++) {
		int m = i;
		double Em = fabs(E[i]);
		for (j = i + 1; j < n; j++) {
			double Ej = fabs(E[j]);
			m = (Em < Ej) ? j : m;
			Em = (Em < Ej) ? Ej : Em;
		}
		if (m != i) {
			int l;
			double b = E[i];
			E[i] = E[m];
			E[m] = b;
			for (j = 0, k = i * n, l = m * n; j < n; j++, k++, l++) {
				b = V[k];
				V[k] = V[l];
				V[l] = b;
			}
		}
	}
	*/
	

	return 0;
}


void matrix_inverse(double T[4][4], double TIN[4][4]) {

	TIN[0][0] = T[0][0];
	TIN[0][1] = T[1][0];
	TIN[0][2] = T[2][0];
	TIN[0][3] = -(T[0][0]*T[0][3] + T[1][0]*T[1][3] + T[2][0]*T[2][3]);

	TIN[1][0] = T[0][1];
	TIN[1][1] = T[1][1];
	TIN[1][2] = T[2][1];
	TIN[1][3] = -(T[0][1]*T[0][3] + T[1][1]*T[1][3] + T[2][1]*T[2][3]);

	TIN[2][0] = T[0][2];
	TIN[2][1] = T[1][2];
	TIN[2][2] = T[2][2];
	TIN[2][3] = -(T[0][2]*T[0][3] + T[1][2]*T[1][3] + T[2][2]*T[2][3]);

	TIN[3][0] = 0;
	TIN[3][1] = 0;
	TIN[3][2] = 0;
	TIN[3][3] = 1;
}


void matrix_mult(double vector[4], double matrix[4][4], double result[4]) {

	result[0] = matrix[0][0] * vector[0] +
				matrix[0][1] * vector[1] +
				matrix[0][2] * vector[2] +
				matrix[0][3] * vector[3];

	result[1] = matrix[1][0] * vector[0] +
				matrix[1][1] * vector[1] +
				matrix[1][2] * vector[2] +
				matrix[1][3] * vector[3];

	result[2] = matrix[2][0] * vector[0] +
				matrix[2][1] * vector[1] +
				matrix[2][2] * vector[2] +
				matrix[2][3] * vector[3];

	result[3] = matrix[3][0] * vector[0] +
				matrix[3][1] * vector[1] +
				matrix[3][2] * vector[2] +
				matrix[3][3] * vector[3];
}


/* End of file */
#endif