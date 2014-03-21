/*
* Software License Agreement (BSD License)
*
*  Object Pose Estimation (OPE) - www.cse.usf.edu/kkduncan/ope
*  Copyright (c) 2013, Kester Duncan
* 
*  Adapted from code by G. Biegelbauer
*/

#include <cstdlib>
#include <cmath>
#include <ctime>
#include <pcl/point_types.h>
#include <Eigen/Core>

#include "SQTypes.h"
#include "Minimization.h"


namespace ope {

Minimization::Minimization() {
	glochisq = 0.;
	i_am_in_trouble = 0;
	
	for (int i = 0; i < MAX_PARAMS; i++) {
		glatry[i] = 0;
		glbeta[i] = 0;
	}

}


Minimization::~Minimization() {
	
}

	

double Minimization::mylog(double x) {
	if (x > 0) { 
		return(log(x));
	} else {
		return(0);
	}    
}


double Minimization::mypow(double x, double y) {
	if (x >= 0) {
		return(pow(x,y));
	} else {
		return(pow(fabs(x),y));
	}    
}


double Minimization::sqr(double x) {
	return (x * x);
}


double Minimization::errorFunc(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams) {
	pcl::PointXYZ _pnt;
	Eigen::MatrixXd m(3, 3), _m(3, 3);

	m(0, 0) = cos(sqParams.phi) * cos(sqParams.theta) * cos(sqParams.psi) - sin(sqParams.phi) * sin(sqParams.psi);
	m(0, 1) = -cos(sqParams.phi) * cos(sqParams.theta) * sin(sqParams.psi) - sin(sqParams.phi) * cos(sqParams.psi);
	m(0, 2) = cos(sqParams.phi) * sin(sqParams.theta);
	m(1, 0) = sin(sqParams.phi) * cos(sqParams.theta) * cos(sqParams.psi) + cos(sqParams.phi) * sin(sqParams.psi);
	m(1, 1) = -sin(sqParams.phi) * cos(sqParams.theta) * sin(sqParams.psi) + cos(sqParams.phi) * cos(sqParams.psi);
	m(1, 2) = sin(sqParams.phi) * sin(sqParams.theta);
	m(2, 0) = -sin(sqParams.theta) * cos(sqParams.psi);
	m(2, 1) = sin(sqParams.theta) * sin(sqParams.psi);
	m(2, 2) = cos(sqParams.theta);

	_m = m.inverse();

	double f, a, b, c, d, e, error = 0;

	for (int i = 0; i < (int) cloud.points.size(); i++) {
		Eigen::MatrixXd pntMove(3, 1);
		Eigen::MatrixXd pnt(3, 1);
		Eigen::MatrixXd _pnt(3, 1);

		pntMove(0, 0) = (double) cloud.points[i].x - sqParams.px; 
		pntMove(1, 0) = (double) cloud.points[i].y - sqParams.py; 
		pntMove(2, 0) = (double) cloud.points[i].z - sqParams.pz;

		_pnt = _m * pntMove;

		_pnt(0, 0) = _pnt(0, 0) / (sqParams.kx * _pnt(2, 0) / sqParams.a3 + 1);
		_pnt(1, 0) = _pnt(1, 0) / (sqParams.ky * _pnt(2, 0) / sqParams.a3 + 1);

		a = pow(fabs(_pnt(0, 0) / sqParams.a1), 2.0 / sqParams.e2);
		b = pow(fabs(_pnt(1, 0) / sqParams.a2), 2.0 / sqParams.e2);
		c = pow(fabs(_pnt(2, 0) / sqParams.a3), 2.0 / sqParams.e1);
		d = pow(a + b, (double) (sqParams.e2 / sqParams.e1));
		e = d + c;
		f = pow(e, (double) (sqParams.e1));

		double tempError = pow((f - 1), 2.0);
		error += exp(-tempError * 10000);
	}

	error = error / static_cast<double>(cloud.points.size());
	error = sqrt(error);

	return (error);
}


double Minimization::funcs(double x, double y, double z, 
	double rotx[11], double roty[11], double rotz[11], 
	glnparam a, glnparam dFda) {

		double coefx, coefy, coefz, A, B, D, G, F, DA, DB, GD, Gz, R, M, N, nx, ny, nz, nz_p, N1, N2, SAB, CAB;
		glnparam dx, dy, dz, dxb, dyb, dzb, dA, dB, dD, dG, dR;
		double xb, yb, zb, koren;
		double O, Omin, Omax, zhat;

		coefx = (rotx[0]*x + roty[0]*y + rotz[0]*z - (a[8]*rotx[0] + a[9]*roty[0] + a[10]*rotz[0]));

		for (int i = 0; i < 5; i++) dx[i] = 0;

		dx[5] = (rotx[3] * x + roty[3] * y + rotz[3] * z - (a[8]*rotx[3] + a[9]*roty[3] + a[10]*rotz[3]));
		dx[6] = (rotx[6] * x + roty[6] * y + rotz[6] * z - (a[8]*rotx[6] + a[9]*roty[6] + a[10]*rotz[6]));
		dx[7] = (rotx[9] * x + roty[9] * y + rotz[9] * z - (a[8]*rotx[9] + a[9]*roty[9] + a[10]*rotz[9]));

		dx[8] = -rotx[0];
		dx[9] = -roty[0];
		dx[10] = -rotz[0];

		for (int i = 11; i < 18; i++) dx[i] = 0;

		coefy = (rotx[1]*x + roty[1]*y + rotz[1]*z - (a[8]*rotx[1] + a[9]*roty[1] + a[10]*rotz[1]));

		for (int i = 0; i < 5; i++) dy[i] = 0;

		dy[5] = (rotx[4] * x + roty[4] * y + rotz[4] * z - (a[8]*rotx[4] + a[9]*roty[4] + a[10]*rotz[4]));
		dy[6] = (rotx[7] * x + roty[7] * y + rotz[7] * z - (a[8]*rotx[7] + a[9]*roty[7] + a[10]*rotz[7]));
		dy[7] = (rotx[10] * x + roty[10] * y + rotz[10] * z - (a[8]*rotx[10] + a[9]*roty[10] + a[10]*rotz[10]));
		dy[8] = -rotx[1];
		dy[9] = -roty[1];
		dy[10] = -rotz[1];

		for (int i = 11; i < 18; i++) dy[i] = 0;

		coefz = (rotx[2]*x + roty[2]*y + rotz[2]*z - (a[8]*rotx[2] + a[9]*roty[2] + a[10]*rotz[2]));

		for (int i = 0; i < 5; i++) dz[i] = 0;

		dz[5] = (rotx[5] * x + roty[5] * y + rotz[5] * z - (a[8]*rotx[5] + a[9]*roty[5] + a[10]*rotz[5]));
		dz[6] = (rotx[8] * x + roty[8] * y + rotz[8] * z - (a[8]*rotx[8] + a[9]*roty[8] + a[10]*rotz[8]));
		dz[7] = 0;
		dz[8] = -rotx[2];
		dz[9] = -roty[2];
		dz[10] = -rotz[2];

		for (int i = 11; i < 18; i++) dz[i] = 0;

		koren = sqrt(sqr(coefx) + sqr(coefy));
		SAB = sin(a[17] - atan2(coefy, coefx));
		CAB = cos(a[17] - atan2(coefy, coefx));

		R = 1/(koren * CAB);

		for (int i = 0; i < 5; i++) 
			dR[i] = 0;

		for (int i = 5; i < 11; i++) 
			dR[i] = (-(coefx * dx[i] + coefy * dy[i]) / CAB - (SAB/sqr(CAB)) * (coefx * dy[i] - coefy * dx[i])) / (koren * sqr(koren));
		
		for (int i = 11; i < 16; i++) 
			dR[i] = 0;

		dR[17] = SAB/(koren * sqr(CAB));


		/**********************************************/

		koren = sqrt(sqr(coefz) + sqr(1/a[16] - 1/R));

		xb = coefx - (1/R - 1/a[16] + koren) * cos(a[17]);

		for (int i = 0; i < 5; i++) 
			dxb[i] = 0;
		
		for (int i = 0; i < 11; i++) 
			dxb[i] = dx[i] - (-dR[i]/sqr(R) + (1.0/koren) * (coefz * dz[i] + (1/a[16] - 1/R) * (dR[i]/sqr(R)))) * cos(a[17]);
		
		for (int i = 11; i < 16; i++) 
			dxb[i] = 0;

		dxb[16] = cos(a[17]) * (-1/sqr(a[16]) + (1/a[16] - 1/R)/(koren * sqr(a[16])));
		dxb[17] = sin(a[17]) * (1/R - 1/a[16] + koren) + cos(a[17]) * dR[17] * (1/sqr(R) - (1/a[16] - 1/R)/(koren*sqr(R)));


		/************************************************/

		yb = coefy - (1/R - 1/a[16] + koren) * sin(a[17]);

		for (int i = 0; i < 5; i++) 
			dyb[i] = 0;

		for (int i = 0; i < 11; i++) 
			dyb[i] = dy[i] - (-dR[i]/sqr(R) + (1.0/koren) * (coefz * dz[i] + (1/a[16] - 1/R) * dR[i]/sqr(R))) * sin(a[17]);
		
		for (int i = 11; i < 16; i++) 
			dyb[i] = 0;

		dyb[16] = sin(a[17]) * (-1/sqr(a[16]) + (1/a[16] - 1/R) / (koren * sqr(a[16])));
		dyb[17] = -cos(a[17]) * (1/R - 1/a[16] + koren) - sin(a[17]) * dR[17] * (-1/sqr(R) + (1/a[16] - 1/R)/(koren* sqr(R)));


		/***************************************************************/

		O = atan2(coefz, 1/a[16] - 1/R);

		zb = O/a[16];

		for (int i = 0; i < 5; i++) 
			dzb[i] = 0;
		
		for (int i = 5; i < 11; i++) 
			dzb[i] = (1/(a[16]*(sqr(1/a[16] - 1/R) + sqr(coefz)))) * ((1/a[16] - 1/R) * dz[i] - coefz * dR[i]/sqr(R));
		
		for (int i = 11; i < 16; i++) 
			dzb[i] = 0;

		dzb[16] = -O/sqr(a[16]) + (1/(a[16] * sqr(a[16]))) * coefz/(sqr(coefz) + sqr(1/a[16] - 1/R));
		dzb[17] = 0;

		/******************************************************************/

		A = xb/(a[0]*(a[11]*zb/a[2] + 1));

		dA[0] = -A/a[0];
		dA[1] = 0;
		dA[2] = (xb/a[0]) * (1.0/sqr(a[11]*zb/a[2] + 1))*a[11]*zb/sqr(a[2]);

		for(int i = 3; i < 11; i++) 
			dA[i] = dxb[i]/(a[0]*(a[11]*zb/a[2] + 1)) - xb * a[11] * dzb[i]/(a[0]*sqr(a[11] * zb/a[2] + 1)*a[2]);

		dA[11] = -xb * zb/(a[0]*sqr(a[11]*zb/a[2] + 1)*a[2]);
		dA[12] = 0;

		for(int i = 13; i < 18; i++) 
			dA[i] = dxb[i]/(a[0]*(a[11]*zb/a[2] + 1)) - xb * a[11] * dzb[i]/(a[0]*sqr(a[11] * zb/a[2] + 1)*a[2]);

		B = yb/(a[1]*(a[12]*zb/a[2] + 1));

		dB[0] = 0;
		dB[1] = -B/a[1];
		dB[2] = (yb/a[1]) * (1.0/sqr(a[12]*zb/a[2] + 1))*a[12]*zb/sqr(a[2]);

		for(int i = 3; i < 11; i++) 
			dB[i] = dyb[i]/(a[1]*(a[12]*zb/a[2] + 1)) - yb * a[12] * dzb[i]/(a[1]*sqr(a[12] * zb/a[2] + 1)*a[2]);

		dB[11] = 0;
		dB[12] = -yb * zb/(a[1]*sqr(a[12]*zb/a[2] + 1)*a[2]);

		for(int i = 13; i < 18; i++) 
			dB[i] = dyb[i]/(a[1]*(a[12]*zb/a[2] + 1)) - yb * a[12] * dzb[i]/(a[1]*sqr(a[12] * zb/a[2] + 1)*a[2]);


		/*********************************************/
		DA = mypow(sqr(A), 1.0/a[4]);   DB = mypow(sqr(B), 1.0/a[4]);
		D = DA + DB;

		for(int i = 0; i < 3; i++)
			dD[i] = (2/a[4]) * (DA/A) * dA[i] + (2/a[4]) * (DB/B) * dB[i];

		dD[3] = 0;
		dD[4] = -DA * mylog(sqr(A))/sqr(a[4]) - DB * mylog(sqr(B))/sqr(a[4]);

		for(int i = 5; i < 18; i++)
			dD[i] = (2/a[4]) * (DA/A) * dA[i] + (2/a[4]) * (DB/B) * dB[i];

		/**********************************************/
		GD = mypow(D, a[4]/a[3]);   Gz = mypow(sqr(zb/a[2]), 1.0/a[3]);
		G = GD + Gz;

		for(int i = 0; i < 2; i++)
			dG[i] = (a[4]/a[3]) * (GD/D) * dD[i];

		dG[2] = (a[4]/a[3]) * (GD/D) * dD[2] - (2/a[3]) * (a[2]/zb)* Gz * (zb/sqr(a[2]));
		dG[3] = GD * (-a[4] * (1/sqr(a[3])) * mylog(D) + (a[4]/a[3])*(1/D) * dD[3]) - (1.0/sqr(a[3])) * Gz * mylog(sqr(zb/a[2]));
		dG[4] = GD * ((1/a[3]) * mylog(D) + (a[4]/a[3])*(1/D) * dD[4]);

		for(int i = 5; i < 18; i++)
			dG[i] = (a[4]/a[3]) * (GD/D) * dD[i] + (2/a[3]) * (Gz/zb) * dzb[i];

		/***********************************************/
		F = mypow(G, a[3]);
		for(int i= 0; i < 3; i++)
			dFda[i] = a[3] * (F/G) * dG[i];

		dFda[3] = F * (mylog(G) + (a[3]/G) * dG[3]);

		for(int i = 4; i < 18; i++)
			dFda[i] = a[3] * (F/G) * dG[i];

		F = F - 1;

		dFda[0] = sqrt(a[1] * a[2]) * (F/(2*sqrt(a[0])) + sqrt(a[0]) * dFda[0]);
		dFda[1] = sqrt(a[0] * a[2]) * (F/(2*sqrt(a[1])) + sqrt(a[1]) * dFda[1]);
		dFda[2] = sqrt(a[0] * a[1]) * (F/(2*sqrt(a[2])) + sqrt(a[2]) * dFda[2]);

		for(int i = 3; i < 18; i++) 
			dFda[i] = dFda[i] * sqrt(a[0] * a[1] * a[2]);

		F = sqrt(a[0] * a[1] * a[2]) * F;

		if(F < 0) {
			return(1 * F);
		} else {
			return(F);
		}
}


void Minimization::precomp_rot(glnparam a, double rotx[11], double roty[11], double rotz[11]) {
	rotx[0] = (cos(a[5])*cos(a[6])*cos(a[7]) - sin(a[5])*sin(a[7]));
	rotx[1] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
	rotx[2] = (cos(a[5])*sin(a[6]));

	rotx[3] = (-sin(a[5])*cos(a[6])*cos(a[7]) - cos(a[5])*sin(a[7]));
	rotx[4] = (sin(a[5])*cos(a[6])*sin(a[7]) - cos(a[5])*cos(a[7]));
	rotx[5] = (-sin(a[5])*sin(a[6]));

	rotx[6] = (-cos(a[5])*sin(a[6])*cos(a[7]));
	rotx[7] = (cos(a[5])*sin(a[6])*sin(a[7]));
	rotx[8] = (cos(a[5])*cos(a[6]));

	rotx[9] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
	rotx[10]= (-cos(a[5])*cos(a[6])*cos(a[7]) + sin(a[5])*sin(a[7]));

	roty[0] = (sin(a[5])*cos(a[6])*cos(a[7]) + cos(a[5])*sin(a[7]));
	roty[1] = (-sin(a[5])*cos(a[6])*sin(a[7]) + cos(a[5])*cos(a[7]));
	roty[2] = (sin(a[5])*sin(a[6]));

	roty[3] = (cos(a[5])*cos(a[6])*cos(a[7]) - sin(a[5])*sin(a[7]));
	roty[4] = (-cos(a[5])*cos(a[6])*sin(a[7]) - sin(a[5])*cos(a[7]));
	roty[5] = (cos(a[5])*sin(a[6]));

	roty[6] = (-sin(a[5])*sin(a[6])*cos(a[7]));
	roty[7] = (sin(a[5])*sin(a[6])*sin(a[7]));
	roty[8] = (sin(a[5])*cos(a[6]));

	roty[9] = (-sin(a[5])*cos(a[6])*sin(a[7]) + cos(a[5])*cos(a[7]));
	roty[10]= (-sin(a[5])*cos(a[6])*cos(a[7]) - cos(a[5])*sin(a[7]));

	rotz[0] = (-sin(a[6])*cos(a[7]));
	rotz[1] = (sin(a[6])*sin(a[7]));
	rotz[2] = (cos(a[6]));

	rotz[3] = 0;
	rotz[4] = 0;
	rotz[5] = 0;

	rotz[6] = (-cos(a[6])*cos(a[7]));
	rotz[7] = (cos(a[6])*sin(a[7]));
	rotz[8] = (-sin(a[6]));

	rotz[9] = (sin(a[6])*sin(a[7]));
	rotz[10] = (sin(a[6])*cos(a[7]));

}


double Minimization::mrqcof(const glndata2& x, glndata F, glndata sig, 
	int ndata, glmma a, gllista lista, int mfit, 
	glcovar alpha, glmma beta, int *n_model, 
	int *n_model_acc, double addnoise) {

	double Fmod, wt, sig2i, dF, chisq, threshold;
	glmma dFda;
	double rotx[11], roty[11], rotz[11];

	for(int j = 0; j < mfit; j++) {
		for(int k = 0; k <= j; k++) {
			alpha[j][k] = 0.0;
		}
		beta[j] = 0.0;
	}

	sig2i = 1.0 / (sig[1] * sig[1]);
	precomp_rot(a, rotx, roty, rotz);	
	threshold = 10e30; // 1.068647458e14

	do { 
		chisq = 0.0;
		*n_model = 0;
		for(int i = 0; i < ndata; i++) {
			Fmod = funcs(x[i][0], x[i][1], x[i][2], rotx, roty, rotz, a, dFda);
			dF = F[i] - Fmod;

			if (Fmod < threshold) {
				(*n_model) +=  1;

				for(int j = 0; j < mfit; j++) {
					wt = dFda[lista[j]]*sig2i;

					for(int k = 0; k <= j; k++) {
						alpha[j][k] = alpha[j][k] + wt * dFda[lista[k]];
					}
					beta[j] = beta[j] + dF*wt;
				}
				chisq = chisq + dF*dF*sig2i;
			}			

			/*****  The iteration cannot be accepted  *********/

			if(chisq / (*n_model_acc) > addnoise) break;
		}
		threshold = threshold * 2;

		if(chisq / (*n_model_acc) > addnoise) break;

	} while ((double)(*n_model) / (*n_model_acc) < 0.95);

	for(int j = 0; j < mfit; j++) {
		for(int k = 0; k <= j ; k++) {
			alpha[j][k] = alpha[j][k] * chisq;
		}
		beta[j] = beta[j] * chisq;
	}

	chisq = chisq / (*n_model);

	for(int j = 1; j < mfit; j++) {
		for(int k = 0; k <= j-1; k++) {
			alpha[k][j] = alpha[j][k]; 
		}
	}

	return(chisq);
}


int Minimization::gaussj(glcovar a, int n, glcovar b) {
	double big, dum, pivinv;
	int icol, irow, m;

	glnp indxc, indxr, ipiv;

	m = 1;

	for(int j = 0; j < n; j++) ipiv[j] = 0;

	for(int i = 0; i < n; i++) {
		big = 0.0;
		for(int j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for(int k = 0; k < n ; k++) {
					if(ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					}
					else if (ipiv[k] > 1) {
						//printf("pause 1 in GAUSSJ - singular matrix\n");
						return(0);			 
					}
				}
			}
		}

		ipiv[icol] = ipiv[icol] + 1;
		if (irow != icol) { 
			for(int l = 0; l < n; l++) {
				dum = a[irow][l];
				a[irow][l] = a[icol][l];
				a[icol][l] = dum;
			}
			for(int l = 0; l < m; l++) {
				dum = b[irow][l];
				b[irow][l] = b[icol][l];
				b[icol][l] = dum;
			}
		}

		indxr[i] = irow;
		indxc[i] = icol;
		//if (a[icol][icol] == 0.0) printf("pause 2 in GAUSSJ - singular matrix\n");

		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for(int l = 0; l < n; l++) a[icol][l] = a[icol][l] *pivinv;
		for(int l = 0; l < m; l ++) b[icol][l] = b[icol][l] * pivinv;
		for(int ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for(int l = 0; l < n; l++) a[ll][l] = a[ll][l] - a[icol][l]* dum;
				for(int l = 0; l < m; l++) b[ll][l] = b[ll][l] - b[icol][l]*dum;
			}
		}
	}

	for(int l = n-1; l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for(int k = 0; k < n; k++) {
				dum = a[k][indxr[l]];
				a[k][indxr[l]] = a[k][indxc[l]];
				a[k][indxc[l]] = dum;
			}
		}
	}

	return(1);
}


double Minimization::mrqmin_init(const glndata2& x, glndata F, glndata sig, int ndata, 
	glmma a, gllista lista, int mfit, glcovar alpha, int nca, int *n_model_acc) { 
	
	int kk, ihit, n_model;
	double alamda, addnoise;

	kk = mfit;
	for(int j = 0; j < nca; j++) {
		ihit = 0;
		for(int k = 0; k < mfit; k++) {
			if (lista[k] == j) ihit = ihit + 1;
		}

		if (ihit == 0) {
			lista[kk] = j; 
			kk = kk + 1;
		}
		else if (ihit > 1) {
			printf("WARNING: pause 1 in routine MRQMIN\n");
			printf("WARNING: Improper permutation in LISTA\n");
		}
	}

	if(kk != nca) {
		printf("WARNING: pause 2 in routine MRQMIN\n");
		printf("WARNING:Improper permutation in LISTA\n");
	}

	alamda = 1;
	addnoise = 1e20; // 485165195.4
	glochisq = mrqcof(x, F, sig, ndata, a, lista, mfit, alpha, glbeta, &n_model, n_model_acc, addnoise);
	*n_model_acc = n_model;

	return(alamda);
}


double Minimization::mrqmin(SQParameters& prm, const glndata2& x, glndata F, glndata sig, 
	int ndata, glmma a, gllista lista, int mfit, glcovar covar, glcovar alpha, 
	double alamda, int *n_model_acc) {
	
	srand(time(NULL));
	double chisq;
	double poisson, addnoise;
	glmma da;
	glcovar oneda;
	int n_model;

	for(int j = 0; j < mfit; j++) {
		for(int k = 0; k < mfit; k++) covar[j][k] = alpha[j][k];
		covar[j][j] = alpha[j][j] * (0.1 + alamda);
		oneda[j][0] = glbeta[j];
	}

	if (gaussj(covar, mfit, oneda) == 0) {
		i_am_in_trouble = 1;
		return(0);
	}

	for(int j = 0; j < mfit; j++) {
		da[j] = oneda[j][0];
	}

	for (int j = 0; j < mfit; j++) {
		if (lista[j] == 0) {
			if (prm.min.a1.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.a1.lowerBound) {
					glatry[lista[j]] = prm.min.a1.lowerBound;
				} else if(a[lista[j]] + da[j] >= prm.min.a1.upperBound) {
					glatry[lista[j]] = prm.min.a1.upperBound;
				} else {
					glatry[lista[j]] = a[lista[j]] + da[j];
				}
			}
			else {
				if (a[lista[j]] + da[j] < 1.0) glatry[lista[j]] = 1.0;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 1) {
			if (prm.min.a2.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.a2.lowerBound) glatry[lista[j]] = prm.min.a2.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.a2.upperBound) glatry[lista[j]] = prm.min.a2.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if (a[lista[j]] + da[j] < 1.0) glatry[lista[j]] = 1.0;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 2) {
			if (prm.min.a3.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.a3.lowerBound) glatry[lista[j]] = prm.min.a3.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.a3.upperBound) glatry[lista[j]] = prm.min.a3.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if (a[lista[j]] + da[j] < 1.0) glatry[lista[j]] = 1.0;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 3) {
			if (prm.min.e1.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.e1.lowerBound) glatry[lista[j]] = prm.min.e1.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.e1.upperBound) glatry[lista[j]] = prm.min.e1.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if(a[lista[j]] + da[j] < 0.1) glatry[lista[j]] = 0.1;
				else if(a[lista[j]] + da[j] > 1.0) glatry[lista[j]] = 1.0;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 4) {
			if (prm.min.e2.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.e2.lowerBound) glatry[lista[j]] = prm.min.e2.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.e2.upperBound) glatry[lista[j]] = prm.min.e2.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if(a[lista[j]] + da[j] < 0.1) glatry[lista[j]] = 0.1;
				else if(a[lista[j]] + da[j] > 1.0) glatry[lista[j]] = 1.0;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 5) {
			if (prm.min.phi.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.phi.lowerBound) glatry[lista[j]] = prm.min.phi.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.phi.upperBound) glatry[lista[j]] = prm.min.phi.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 6) {
			if (prm.min.theta.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.theta.lowerBound) glatry[lista[j]] = prm.min.theta.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.phi.upperBound) glatry[lista[j]] = prm.min.theta.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 7) {
			if (prm.min.psi.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.psi.lowerBound) glatry[lista[j]] = prm.min.psi.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.psi.upperBound) glatry[lista[j]] = prm.min.psi.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 8) {
			if (prm.min.px.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.px.lowerBound) glatry[lista[j]] = prm.min.px.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.px.upperBound) glatry[lista[j]] = prm.min.px.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 9) {
			if (prm.min.py.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.py.lowerBound) glatry[lista[j]] = prm.min.py.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.py.upperBound) glatry[lista[j]] = prm.min.py.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 10) {
			if (prm.min.pz.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.pz.lowerBound) glatry[lista[j]] = prm.min.pz.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.pz.upperBound) glatry[lista[j]] = prm.min.pz.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if (lista[j] == 11) {
			if (prm.min.kx.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.kx.lowerBound) glatry[lista[j]] = prm.min.kx.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.kx.upperBound) glatry[lista[j]] = prm.min.kx.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if(a[lista[j]] + da[j] > 1) glatry[lista[j]] = 1;
				else if(a[lista[j]] + da[j] < -1) glatry[lista[j]] = -1;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if (lista[j] == 12) {
			if (prm.min.ky.type == BOUNDED) {
				if(a[lista[j]] + da[j] <= prm.min.ky.lowerBound) glatry[lista[j]] = prm.min.ky.lowerBound;
				else if(a[lista[j]] + da[j] >= prm.min.ky.upperBound) glatry[lista[j]] = prm.min.ky.upperBound;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
			else {
				if(a[lista[j]] + da[j] > 1) glatry[lista[j]] = 1;
				else if(a[lista[j]] + da[j] < -1) glatry[lista[j]] = -1;
				else glatry[lista[j]] = a[lista[j]] + da[j];
			}
		}

		else if(lista[j] == 13) {
			if(a[lista[j]] + da[j] < a[14] || a[lista[j]] + da[j] > a[15]) glatry[lista[j]] = a[lista[j]];
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if(lista[j] == 14) {
			if(a[lista[j]] + da[j] > a[13]) glatry[lista[j]] = a[lista[j]];
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if(lista[j] == 15) {
			if(a[lista[j]] + da[j] < a[13]) glatry[lista[j]] = a[lista[j]];
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if(lista[j] == 16) {
			if(fabs(a[lista[j]] + da[j]) > 1e9) glatry[lista[j]] = a[lista[j]];
			else if(a[lista[j]] + da[j] < 0.0) glatry[lista[j]] = a[lista[j]];
			else glatry[lista[j]] = a[lista[j]] + da[j];
		}

		else if(lista[j] == 17) {
			glatry[lista[j]] = a[lista[j]] + da[j];
			glatry[lista[j]] = glatry[lista[j]] - 2 * PI * (int)(glatry[lista[j]]/(2 * PI));
		}

		else {
			glatry[lista[j]] = a[lista[j]]+ da[j];
		}

	}

	for (int j = mfit; j < MAX_PARAMS; j++) {
		glatry[lista[j]] = a[lista[j]];
	}

	/// Adding Noise
	poisson = ((double) rand()) / RAND_MAX;
	addnoise = glochisq + (glochisq) * fabs(poisson) / 2;

	chisq = mrqcof(x, F, sig, ndata, glatry, lista, mfit, covar, da, &n_model, n_model_acc, addnoise);

	if (chisq < addnoise) {
		if( chisq < addnoise) alamda = alamda / NU;
		glochisq = chisq;
		*n_model_acc = n_model;
		for(int j = 0; j < mfit; j++) {
			glbeta[j] = da[j];
			a[lista[j]] = glatry[lista[j]];
			for(int k = 0; k < mfit; k++) alpha[j][k] = covar[j][k];
		}
	}
	else {
		alamda = NU * alamda;
		chisq = glochisq;
	}

	return(alamda);

}


} /* ope */
