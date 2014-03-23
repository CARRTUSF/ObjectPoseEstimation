/*
* Software License Agreement (BSD License)
*
*  Object Pose Estimation (OPE) - www.cse.usf.edu/kkduncan/ope
*  Copyright (c) 2013, Kester Duncan
*
*  All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*
*	\file	Minimization.h
*	\brief	Defines functions for Levenberg-Marquadt minimization etc.
*	\author	G. Biegelbauer, restructured by Kester Duncan
*/
#pragma once
#ifndef MINIMIZATION_H__
#define MINIMIZATION_H__

//#include <cstdlib>


////////////// Forward declarations to avoid large headers ////////////////
namespace pcl {
	struct PointXYZ;
	template<class T> class PointCloud;

} // pcl

class SQParameters;

//////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////
// Types and variables used by the Levenberg Marquadt minimization procedure
////////////////////////////////////////////////////////////////////////////////////

/// Maximum number of points
#define MAX_POINTS 180000

/// Maximum number of parameters to minimize
#define MAX_PARAMS 18

/// Parameter used during the minimization process
#define	NU 10.0

/// Literal value of PI
#define PI 3.141592654
	



namespace ope {

	
/**
 * \brief Defines a 3D point that would be used by the minimization routines herein
 */
struct MPoint {
public:
	double point[3];

	MPoint() {
		point[0] = 0.0;
		point[1] = 0.0;
		point[2] = 0.0;
	}

	inline double& operator[] (const int& idx) {
		return point[idx];
	}

	inline double operator[] (const int& idx) const {
		return point[idx];
	}

};


typedef double				glmma[MAX_PARAMS];
typedef glmma				glnparam;
typedef int					gllista[MAX_PARAMS];
typedef double				glcovar[MAX_PARAMS][MAX_PARAMS];
typedef glcovar				glnalbynal;
typedef glcovar				glncabynca;
typedef glcovar				glnpbynp;
typedef glcovar				glnpbymp;
typedef std::vector<double>	glndata;
typedef std::vector<MPoint>	glndata2;
typedef gllista				glnp;
typedef double				glmartrix[4][4];


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////





/**
 * \brief Defines the minimization procedures
 */
struct Minimization {
	/**
	 * Properties used throughout pose estimation process
	 */
	double	glochisq;
	glmma	glatry;
	glmma	glbeta;
	int		i_am_in_trouble;


	Minimization();
	~Minimization();


	/**
	 * \brief Calculates the base-10 logarithm of x
	 */
	double mylog(double x);


	/**
	 * \brief Raises x to the power of y
	 */
	double mypow(double x, double y);


	/**
	 * \brief Calculates the square root of x
	 */
	double sqr(double x);


	/**
	 * \brief Calculates the average error of fit of superquadric parameters.
	 * \details The error of fit is the mean of the algebraic distance calculated by the inside-outside function.
	 */
	double errorFunc(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams);


	/**
	 * \brief Levenberg-Marquadt cost function for superquadric fitting
	 */
	double funcs(double x, double y, double z, 
		double rotx[11], double roty[11], double rotz[11], 
		glnparam a, glnparam dFda);
	
	
	/**
	 * \brief Utility function for computing rotations
	 */
	void precomp_rot(glnparam a, double rotx[11], double roty[11], double rotz[11]);


	/**
	 * \brief Utility function for Levenberg-Marquadt minimization
	 */
	double mrqcof(const glndata2& x, glndata F, glndata sig, int ndata, 
		glmma a, gllista lista, int mfit, glcovar alpha, glmma beta, 
		int *n_model, int *n_model_acc, double addnoise);


	/**
	 * \brief Utility function for Levenberg-Marquadt minimization
	 */
	int gaussj(glcovar a, int n, glcovar b);


	/**
	 * \brief Initialize the Levenberg-Marquadt minimization
	 */
	double mrqmin_init(const glndata2& x, glndata F, glndata sig, int ndata, 
		glmma a, gllista lista, int mfit, glcovar alpha, int nca, int *n_model_acc);


	/**
	 * \brief Levenberg Marquadt minimization main procedure
	 */
	double mrqmin(SQParameters& prm, const glndata2& x, glndata F, glndata sig, 
		int ndata, glmma a, gllista lista, int mfit, glcovar covar, glcovar alpha, 
		double alamda, int *n_model_acc);


};


} /* ope */

#endif /* MINIMIZATION_H__ */