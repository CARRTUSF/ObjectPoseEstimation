/*
 * Software License Agreement (BSD License)
 *
 *  Object Pose Estimation (OPE) - www.cse.usf.edu/kkduncan/ope
 *  Copyright (c) 2013, Kester Duncan
 *
 */

#include <pcl/point_types.h>
#include <Eigen/Core>
#include "EigenSolver.h"
#include "OPESettings.h"
#include "SQTypes.h"
#include "SQFitting.h"


#ifndef M_PI
#define M_PI 3.141592654
#endif


namespace ope {

SQFitting::SQFitting() {
	
}


SQFitting::~SQFitting() {
	
}


void SQFitting::initializeMinimizationParameters (SQParameters& sqParams) {
	/// Assign size parameter values based on their bounds properties
	if (sqParams.min.a1.type == UNCHANGED) {
		sqParams.a1 = sqParams.min.a1.value;
	}
	if (sqParams.min.a2.type == UNCHANGED) {
		sqParams.a2 = sqParams.min.a2.value;
	}
	if (sqParams.min.a3.type == UNCHANGED) {
		sqParams.a3 = sqParams.min.a3.value;
	}

	if (sqParams.min.a1.type == BOUNDED) {
		sqParams.a1 = (sqParams.min.a1.upperBound + sqParams.min.a1.lowerBound) / 2;
	}
	if (sqParams.min.a2.type == BOUNDED) {
		sqParams.a2 = (sqParams.min.a2.upperBound + sqParams.min.a2.lowerBound) / 2;
	}
	if (sqParams.min.a3.type == BOUNDED) {
		sqParams.a3 = (sqParams.min.a3.upperBound + sqParams.min.a3.lowerBound) / 2;
	}

	if (sqParams.min.a1.type == NOT_USED || sqParams.min.a1.type == UNLIMITED) {
		sqParams.a1 = 1;
	}
	if (sqParams.min.a2.type == NOT_USED || sqParams.min.a2.type == UNLIMITED) {
		sqParams.a2 = 1;
	}
	if (sqParams.min.a3.type == NOT_USED || sqParams.min.a3.type == UNLIMITED) {
		sqParams.a3 = 1;
	}	
}



int SQFitting::estimateParameters(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams, const int& eigenVector) {
	double t[4][4] = {0};

	// Check if cloud contains enough data points
	int len = (int) cloud.points.size();
	if (len < 13) { // changed from 13
		return -1;
	}

	// Find center of gravity
	double a, b, c, x, y, z, xx, yy, zz, cxy, cxz, cyz, count = 0.0;
	
	x = y = z = 0; 
	xx = yy = zz = cxy = cyz = cxz = 0; 
	count = 0;

	for(int i = 0; i < len; i++) {
		pcl::PointXYZ pnt = cloud.points[i];
		a = (double) pnt.x;
		b = (double) pnt.y;
		c = (double) pnt.z;

		x += a;
		y += b;
		z += c;

		xx += a*a;
		yy += b*b;
		zz += c*c;

		cxz += a*c;
		cyz += b*c;
		cxy += a*b;

		count++;
	}

	// Compute center of gravity, covariance and variance
	double avgx, avgy, avgz;
	avgx = x/count;
	avgy = y/count;
	avgz = z/count;

	xx  /= count;
	yy  /= count;
	zz  /= count;

	cxz /= count;
	cyz /= count;
	cxy /= count;

    xx  -= (avgx * avgx);
    yy  -= (avgy * avgy);
    zz  -= (avgz * avgz);

    cxy -= (avgx * avgy);
    cxz -= (avgx * avgz);
    cyz -= (avgy * avgz);

	// Build covariance matrix
	double aa[3][3];
	aa[0][0] =  xx;   aa[0][1] = cxy; aa[0][2] = cxz;
	aa[1][0] = cxy;   aa[1][1] =  yy; aa[1][2] = cyz;
	aa[2][0] = cxz;   aa[2][1] = cyz; aa[2][2] =  zz;

	// Compute eigenvectors and eigenvalues
	double	d[3], m_brain[3][3], m_inverse[3][3];
	double lambdas[3];
	double vectors[3][3];
	int    nrot;

	/*
	 * Using the covariance matrix, we get its eigenvalues
	 * and eigenvectors
	 */
	cvJacobiEigens_64d((double*) aa, (double*) m_brain, (double*) d, 3, 0.0);

	// Find the vector with the largest eigenvalue
	double max_eigen = -1000000000.0;
	int vectorIndex;
	
	if(d[0] > max_eigen) {
		max_eigen = d[0];
		vectorIndex = 0;
	}
	if(d[1] > max_eigen) {
		max_eigen = d[1]; 
		vectorIndex = 1;
	}
	if(d[2] > max_eigen) {
		max_eigen = d[2]; 
		vectorIndex = 2;
	}
	
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			m_inverse[j][k] = m_brain[j][k];
		}
	}

	// Aligning max eigenvector with z axis
#if 0

	for (int j = 0; j < 3; j++) { 
		m_inverse[j][2] = m_brain[j][eigenVector];
		m_inverse[j][eigenVector] = m_brain[j][2];
	}

	
	if (vectorIndex != 2) {
		for (int j = 0; j < 3; j++) {
			m_inverse[j][vectorIndex] = -m_inverse[j][vectorIndex];
		}
	}	

#endif

	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			m_brain[j][k] = m_inverse[j][k];
		}
	}
	
	// Build transformation matrix (rotation & translation)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++){
			t[i][j] = m_brain[i][j];
	    }
	}
	
	t[0][3] = avgx;
	t[1][3] = avgy;
	t[2][3] = avgz;
	t[3][3] =  1;

	t[3][0] = t[3][1] = t[3][2] = 0.0;

	// Set / Calculate initial parameter
	sqParams.e1 = 1;
	sqParams.e2 = 1;
	sqParams.kx = 0;
	sqParams.ky = 0;
	
	double xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = ymin = zmin = 100000.00;
	xmax = ymax = zmax = -100000.00;

	double tInv[4][4];
	matrix_inverse(t, tInv);

	double vector[4];
	for(int i = 0; i < len; i++) {
		pcl::PointXYZ pnt = cloud.points[i];
		
		vector[0] = (double) pnt.x;
		vector[1] = (double) pnt.y;
		vector[2] = (double) pnt.z;
		vector[3] = 1;

		double result[4];
		matrix_mult(vector, tInv, result);
		a = result[0]; 
		b = result[1]; 
		c = result[2];

		if(xmin > a)  xmin = a;     
		if(xmax < a)  xmax = a;

		if(ymin > b)  ymin = b;     
		if(ymax < b)  ymax = b;

		if(zmin > c)  zmin = c;     
		if(zmax < c)  zmax = c;
	}

	double r1, r2, r3;

	r1 = atan2(t[1][2], t[0][2]);
	r1 = r1 + M_PI;
	r2 = atan2(cos(r1) * t[0][2] + sin(r1) * t[1][2], t[2][2]);	
	r3 = atan2(-sin(r1) * t[0][0] + cos(r1) * t[1][0], -sin(r1) * t[0][1] + cos(r1) * t[1][1]);
	
	sqParams.phi = r1;
	sqParams.theta = r2;
	sqParams.psi = r3;

	sqParams.px = t[0][3] + (xmax + xmin) / 2;   
	sqParams.py = t[1][3] + (ymax + ymin) / 2;
	sqParams.pz = t[2][3] + (zmax + zmin) / 2;
 
	double ta1, ta2, ta3;

	ta1 = (xmax - xmin)/2;
	ta2 = (ymax - ymin)/2;
	ta3 = (zmax - zmin)/2;

	sqParams.a1 = ta1;
	sqParams.a2 = ta2;
	sqParams.a3 = ta3;	

	return 0;
}



int SQFitting::estimateParametersNew(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams, const int& eigenVector) {
	double t[4][4] = {0};

	// Check if cloud contains enough data points
	int len = (int) cloud.points.size();
	if (len < 13) { // changed from 13
		return -1;
	}

	// Find center of gravity
	double a, b, c, x, y, z, xx, yy, zz, cxy, cxz, cyz, count;
	
	x = y = z = 0; 
	xx = yy = zz = cxy = cyz = cxz = 0; 
	count = 0;

	for(int i = 0; i < len; i++) {
		pcl::PointXYZ pnt = cloud.points[i];
		a = (double) pnt.x;
		b = (double) pnt.y;
		c = (double) pnt.z;

		x += a;
		y += b;
		z += c;

		xx += a*a;
		yy += b*b;
		zz += c*c;

		cxz += a*c;
		cyz += b*c;
		cxy += a*b;

		count++;
	}

	// Compute center of gravity, covariance and variance
	double avgx, avgy, avgz;
	avgx = x/count;
	avgy = y/count;
	avgz = z/count;

	xx  /= count;
	yy  /= count;
	zz  /= count;

	cxz /= count;
	cyz /= count;
	cxy /= count;

    xx  -= (avgx * avgx);
    yy  -= (avgy * avgy);
    zz  -= (avgz * avgz);

    cxy -= (avgx * avgy);
    cxz -= (avgx * avgz);
    cyz -= (avgy * avgz);

	// Build covariance matrix
	double aa[3][3];
	aa[0][0] =  xx;   aa[0][1] = cxy; aa[0][2] = cxz;
	aa[1][0] = cxy;   aa[1][1] =  yy; aa[1][2] = cyz;
	aa[2][0] = cxz;   aa[2][1] = cyz; aa[2][2] =  zz;

	// Compute eigenvectors and eigenvalues
	double	d[3], m_brain[3][3], m_inverse[3][3];
	double lambdas[3], vectors[3][3];
	int     nrot;

	/*
	 * Using the covariance matrix, we get its eigenvalues
	 * and eigenvectors
	 */
	cvJacobiEigens_64d((double*) aa, (double*) m_brain, (double*) d, 3, 0.0);

	// Find the vector with the largest eigenvalue
	double max_eigen = -1000000000.0;
	int vectorIndex;
	int longestAxis;
	int longestAxisLambda;

#if 0
	if (abs(lambdas[2] - lambdas[1]) < abs(lambdas[1] - lambdas[0])) {
		longestAxisLambda = 0;
	} else {
		longestAxisLambda = 2;
	}

	for (int i = 0; i < 3; i++) {
		if (vectors[0][longestAxisLambda] == m_brain[0][i] &&
			vectors[1][longestAxisLambda] == m_brain[1][i] &&
			vectors[2][longestAxisLambda] == m_brain[2][i]) {
				vectorIndex = i;
		}
	}

#endif

	if(d[0] > max_eigen) {
		max_eigen = d[0];
		vectorIndex = 0;
	}
	if(d[1] > max_eigen) {
		max_eigen = d[1]; 
		vectorIndex = 1;
	}
	if(d[2] > max_eigen) {
		max_eigen = d[2]; 
		vectorIndex = 2;
	}
	
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			m_inverse[j][k] = m_brain[j][k];
		}
	}

	// Aligning max eigenvector with z axis
#if 0
	for (int j = 0; j < 3; j++) { 
		m_inverse[j][2] = m_brain[j][eigenVector];
		m_inverse[j][eigenVector] = m_brain[j][2];
	}

	if (vectorIndex != 2) {
		for (int j = 0; j < 3; j++) {
			m_inverse[j][vectorIndex] = -m_inverse[j][vectorIndex];
		}
	}	

	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			m_brain[j][k] = m_inverse[j][k];
		}
	}

#endif
	
	// Build transformation matrix (rotation & translation)
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++){
			t[i][j] = m_brain[i][j];
	    }
	}
	
	t[0][3] = avgx;
	t[1][3] = avgy;
	t[2][3] = avgz;
	t[3][3] =  1;

	t[3][0] = t[3][1] = t[3][2] = 0.0;

	// Set / Calculate initial parameter
	sqParams.e1 = 1;
	sqParams.e2 = 1;
	sqParams.kx = 0;
	sqParams.ky = 0;
	
	double xmin, ymin, zmin, xmax, ymax, zmax;
	xmin = ymin = zmin = 100000.00;
	xmax = ymax = zmax = -100000.00;

	double tInv[4][4];
	matrix_inverse(t, tInv);

	double vector[4];
	for(int i = 0; i < len; i++) {
		pcl::PointXYZ pnt = cloud.points[i];
		
		vector[0] = (double) pnt.x;
		vector[1] = (double) pnt.y;
		vector[2] = (double) pnt.z;
		vector[3] = 1;

		double result[4];
		matrix_mult(vector, tInv, result);
		a = result[0]; 
		b = result[1]; 
		c = result[2];

		if(xmin > a)  xmin = a;     
		if(xmax < a)  xmax = a;

		if(ymin > b)  ymin = b;     
		if(ymax < b)  ymax = b;

		if(zmin > c)  zmin = c;     
		if(zmax < c)  zmax = c;
	}

	double r1, r2, r3;

	r1 = atan2(t[1][2], t[0][2]);
	r1 = r1 + M_PI;
	r2 = atan2(cos(r1) * t[0][2] + sin(r1) * t[1][2], t[2][2]);	
	r3 = atan2(-sin(r1) * t[0][0] + cos(r1) * t[1][0], -sin(r1) * t[0][1] + cos(r1) * t[1][1]);
	
	sqParams.phi = r1;
	sqParams.theta = r2;
	sqParams.psi = r3;

	sqParams.px = t[0][3] + (xmax + xmin) / 2;   
	sqParams.py = t[1][3] + (ymax + ymin) / 2;
	sqParams.pz = t[2][3] + (zmax + zmin) / 2;
 
	double ta1, ta2, ta3;

	ta1 = (xmax - xmin)/2;
	ta2 = (ymax - ymin)/2;
	ta3 = (zmax - zmin)/2;

	sqParams.a1 = ta1;
	sqParams.a2 = ta2;
	sqParams.a3 = ta3;	

	return vectorIndex;

}



void SQFitting::estimateInitialParameters(pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams) {
	double minError = 1000000;
	int maxEigenVector = 0;
	int retVal;
	
	for (int i = 0; i < 3; i++) {
		retVal = estimateParameters(cloud, sqParams, i);
		if (retVal == -1) {
			printf("Error: At least 13 data points must support the fit!\n");
		}
		
		double error = 0;
		error = minProcess.errorFunc(cloud, sqParams);
		
		if (error < minError) {
			minError = error;
			maxEigenVector = i;
		}
	}

	retVal = estimateParameters(cloud, sqParams, maxEigenVector);
}


double SQFitting::qualityOfFit(pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams) {
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

	error = sqrt(error);

	return (error);
}


double SQFitting::recoverParameters(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& prm) {
	int gliset;
	double gloldm;
	
	int npt, mfit, n_model_acc, iter;
	double alamda, old_chisq;

	glmma a;
	glndata F;
	glndata sig;
	glndata2 xw;
	gllista lista;
	glcovar covar, alpha;
	
	gliset = 0;
	gloldm = -1.0;
	minProcess.i_am_in_trouble = 0;
	
	a[0]  = prm.a1;
	a[1]  = prm.a2;
	a[2]  = prm.a3;
	a[3]  = prm.e1;
	a[4]  = prm.e2;
	a[5]  = prm.phi;
	a[6]  = prm.theta;
	a[7]  = prm.psi;
	a[8]  = prm.px;
	a[9]  = prm.py;
	a[10] = prm.pz;
	a[11] = prm.kx;
	a[12] = prm.ky;
	a[13] = 0;           /* center of bending */
	a[14] = -500;        /* z min */
	a[15] = 500;         /* z max */
	a[16] = 0.0000010;   /* k or 1/k ?? */
	a[17] = 0;           /* alpha angle */
	
	for (int i = 0; i < MAX_PARAMS; i++) lista[i] = i;

	mfit = 13;
	int index = 0;

	if (prm.min.a1.type == NOT_USED || prm.min.a1.type == UNCHANGED) mfit--;
	else {
		lista[index] = 0;
		index++;
	}	
	if (prm.min.a1.type == UNCHANGED) a[0] = prm.min.a1.value;
	else if (prm.min.a1.type == BOUNDED) a[0] = (prm.min.a1.upperBound + prm.min.a1.lowerBound) / 2;
	
	if (prm.min.a2.type == NOT_USED || prm.min.a2.type == UNCHANGED) mfit--;
	else {
		lista[index] = 1;
		index++;
	}
	if (prm.min.a2.type == UNCHANGED) a[1] = prm.min.a2.value;
	else if (prm.min.a2.type == BOUNDED) a[1] = (prm.min.a2.upperBound + prm.min.a2.lowerBound) / 2;
	
	if (prm.min.a3.type == NOT_USED || prm.min.a3.type == UNCHANGED) mfit--;
	else {
		lista[index] = 2;
		index++;
	}
	if (prm.min.a3.type == UNCHANGED) a[2] = prm.min.a3.value;
	else if (prm.min.a3.type == BOUNDED) a[2] = (prm.min.a3.upperBound + prm.min.a3.lowerBound) / 2;
	
	if (prm.min.e1.type == NOT_USED || prm.min.e1.type == UNCHANGED) mfit--;
	else {
		lista[index] = 3;
		index++;
	}
	if (prm.min.e1.type == UNCHANGED) a[3] = prm.min.e1.value;
	else if (prm.min.e1.type == BOUNDED) a[3] = (prm.min.e1.upperBound + prm.min.e1.lowerBound) / 2;
	
	if (prm.min.e2.type == NOT_USED || prm.min.e2.type == UNCHANGED) mfit--;
	else {
		lista[index] = 4;
		index++;
	}
	if (prm.min.e2.type == UNCHANGED) a[4] = prm.min.e2.value;
	else if (prm.min.e2.type == BOUNDED) a[4] = (prm.min.e2.upperBound + prm.min.e2.lowerBound) / 2;
	
	if (prm.min.phi.type == NOT_USED || prm.min.phi.type == UNCHANGED) mfit--;
	else {
		lista[index] = 5;
		index++;
	}
	if (prm.min.phi.type == UNCHANGED) a[5] = prm.min.phi.value;
	else if (prm.min.phi.type == BOUNDED) a[5] = (prm.min.phi.upperBound + prm.min.phi.lowerBound) / 2;
	
	if (prm.min.theta.type == NOT_USED || prm.min.theta.type == UNCHANGED) mfit--;
	else {
		lista[index] = 6;
		index++;
	}
	if (prm.min.theta.type == UNCHANGED) a[6] = prm.min.theta.value;
	else if (prm.min.phi.type == BOUNDED) a[6] = (prm.min.theta.upperBound + prm.min.theta.lowerBound) / 2;
	
	if (prm.min.psi.type == NOT_USED || prm.min.psi.type == UNCHANGED) mfit--;
	else {
		lista[index] = 7;
		index++;
	}
	if (prm.min.psi.type == UNCHANGED) a[7] = prm.min.psi.value;
	else if (prm.min.psi.type == BOUNDED) a[7] = (prm.min.psi.upperBound + prm.min.psi.lowerBound) / 2;
	
	if (prm.min.px.type == NOT_USED || prm.min.px.type == UNCHANGED) mfit--;
	else {
		lista[index] = 8;
		index++;
	}
	if (prm.min.px.type == UNCHANGED) a[8] = prm.min.px.value;
	else if (prm.min.px.type == BOUNDED) a[8] = (prm.min.px.upperBound + prm.min.px.lowerBound) / 2;
	
	if (prm.min.py.type == NOT_USED || prm.min.py.type == UNCHANGED) mfit--;
	else {
		lista[index] = 9;
		index++;
	}
	if (prm.min.py.type == UNCHANGED) a[9] = prm.min.py.value;
	else if (prm.min.py.type == BOUNDED) a[9] = (prm.min.py.upperBound + prm.min.py.lowerBound) / 2;
	
	if (prm.min.pz.type == NOT_USED || prm.min.pz.type == UNCHANGED) mfit--;
	else {
		lista[index] = 10;
		index++;
	}
	if (prm.min.pz.type == UNCHANGED) a[10] = prm.min.pz.value;
	else if (prm.min.pz.type == BOUNDED) a[10] = (prm.min.pz.upperBound + prm.min.pz.lowerBound) / 2;
	
	if (prm.min.kx.type == NOT_USED || prm.min.kx.type == UNCHANGED) mfit--;
	else {
		lista[index] = 11;
		index++;
	}
	if (prm.min.kx.type == UNCHANGED) a[11] = prm.min.kx.value;
	else if (prm.min.kx.type == BOUNDED) a[11] = (prm.min.kx.upperBound + prm.min.kx.lowerBound) / 2;
	
	if (prm.min.ky.type == NOT_USED || prm.min.ky.type == UNCHANGED) mfit--;
	else {
		lista[index] = 12;
		index++;
	}
	if (prm.min.ky.type == UNCHANGED) prm.ky = prm.min.ky.value;
	else if (prm.min.ky.type == BOUNDED) a[12] = (prm.min.ky.upperBound + prm.min.ky.lowerBound) / 2;

	// Prepare surface points
	npt = (int) cloud.points.size();

	if (npt > MAX_POINTS) {
		npt = MAX_POINTS;
	}

	for (int i = 0; i < npt; i++) {
		pcl::PointXYZ pt = cloud.points[i];
		ope::MPoint p;
		p[0] = (double) pt.x;
		p[1] = (double) pt.y;
		p[2] = (double) pt.z;
		xw.push_back(p);
		F.push_back(0);
		sig.push_back(1);		
		
	}

	n_model_acc = npt;
	minProcess.glochisq = 10e31;
	old_chisq = 1.0;
	iter = 0;

	for (int lmIter = 0; lmIter < prm.min.iterations; lmIter++) {
		if (lmIter == 0) 
			alamda = minProcess.mrqmin_init(xw, F, sig, npt, a, lista, mfit, alpha, MAX_PARAMS, &n_model_acc);
		else 
			alamda = minProcess.mrqmin(prm, xw, F, sig, npt, a, lista, mfit, covar, alpha, alamda, &n_model_acc);

		if (minProcess.i_am_in_trouble) {
			return(-1);
		}

		if(old_chisq != minProcess.glochisq) { 
			old_chisq = minProcess.glochisq;
			iter = iter + 1;
		}

	}
	
	///////// Set Return Values ////////
	// Size Parameters
	if (prm.min.a1.type != NOT_USED && prm.min.a1.type != UNCHANGED) prm.a1 = a[0];
	else if (prm.min.a1.type == BOUNDED) prm.a1 = a[0];
	else prm.a1 = prm.min.a1.value;
	
	if (prm.min.a2.type != NOT_USED && prm.min.a2.type != UNCHANGED) prm.a2 = a[1];
	else if (prm.min.a2.type == BOUNDED) prm.a2 = a[1];
	else prm.a2 = prm.min.a2.value;
	
	if (prm.min.a3.type != NOT_USED && prm.min.a3.type != UNCHANGED) prm.a3 = a[2];
	else if (prm.min.a3.type == BOUNDED) prm.a3 = a[2];
	else prm.a3 = prm.min.a3.value;


	// Shape Parameters
	if (prm.min.e1.type != NOT_USED && prm.min.e1.type != UNCHANGED) prm.e1 = a[3];
	else if (prm.min.e1.type == BOUNDED) prm.e1 = a[3];
	else prm.e1 = prm.min.e1.value;
	
	if (prm.min.e2.type != NOT_USED && prm.min.e2.type != UNCHANGED) prm.e2 = a[4];
	else if (prm.min.e2.type == BOUNDED) prm.e2 = a[4];
	else prm.e2 = prm.min.e2.value;


	// Rotation Parameters
	if (prm.min.phi.type != NOT_USED && prm.min.phi.type != UNCHANGED) prm.phi = a[5];
	else if (prm.min.phi.type == BOUNDED) prm.phi = a[5];
	else prm.phi = prm.min.phi.value;
	
	if (prm.min.theta.type != NOT_USED && prm.min.theta.type != UNCHANGED) prm.theta = a[6];
	else if (prm.min.theta.type == BOUNDED) prm.theta = a[6];
	else prm.theta = prm.min.theta.value;
	
	if (prm.min.psi.type != NOT_USED && prm.min.psi.type != UNCHANGED) prm.psi = a[7];
	else if (prm.min.psi.type == BOUNDED) prm.psi = a[7];
	else prm.psi = prm.min.psi.value;


	// Position (translation) Parameters
	if (prm.min.px.type != NOT_USED && prm.min.px.type != UNCHANGED) prm.px = a[8];
	else if (prm.min.px.type == BOUNDED) prm.px = a[8];
	else prm.px = prm.min.px.value;
	
	if (prm.min.py.type != NOT_USED && prm.min.py.type != UNCHANGED) prm.py = a[9];
	else if (prm.min.py.type == BOUNDED) prm.py = a[9];
	else prm.py = prm.min.py.value;
	
	if (prm.min.pz.type != NOT_USED && prm.min.pz.type != UNCHANGED) prm.pz = a[10];
	else if (prm.min.pz.type == BOUNDED) prm.pz = a[10];
	else prm.pz = prm.min.pz.value;


	// Taper Parameters
	if (prm.min.kx.type != NOT_USED && prm.min.kx.type != UNCHANGED) prm.kx = a[11];
	else if (prm.min.kx.type == BOUNDED) prm.kx = a[11];
	else prm.kx = prm.min.kx.value;
	
	if (prm.min.ky.type != NOT_USED && prm.min.ky.type != UNCHANGED) prm.ky = a[12];
	else if (prm.min.ky.type == BOUNDED) prm.ky = a[12];
	else prm.ky = prm.min.ky.value;

	double error = 0;
	error = minProcess.errorFunc(cloud, prm);	

	return (error);	
}



void SQFitting::performShapeFitting(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& initParams, SQParameters& bestParams) {
	pcl::PointCloud<pcl::PointXYZ> tempSQCloud;
	SQParameters unrefinedParams, refinedParams;
	double error, minError = 1.0e13;;
	int ret = -1, eigenVect;
	int iterations;

	// Get best eigenvector
	iterations = initParams.min.iterations;
	initParams.min.iterations = 20;
		
	for (int i = 0; i < 3; i++) {
		ret = estimateParameters(cloud, initParams, i);
		if (ret == -1) {
			printf("Error: Cannot perform superquadric fitting!\n");
		}
		
		error = recoverParameters(cloud, initParams);
		if (error < minError) {
			minError = error;
			eigenVect = i;
		}
	}

	initParams.min.iterations = iterations;
	initParams.principalAxis = eigenVect;
	
	ret = estimateParameters(cloud, initParams, eigenVect);
	if (ret == -1) {
		printf("Error: Cannot perform superquadric fitting!\n");
		error = 100000;
	} else {
		error = recoverParameters(cloud, initParams);
	}

	initParams.copyTo(bestParams);
}
	

} // ope


