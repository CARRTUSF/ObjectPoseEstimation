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
 *	\file	SQFitting.h
 *	\brief	Performs superquadric fitting on a point cloud
 *	\author	Kester Duncan
 *	\note Adapted from G. Biegelbauer's implementation
 *	\details Defines all the functions necessary for Superquadric fitting
 */
#pragma once
#ifndef SQ_FITTING_H__
#define SQ_FITTING_H__


#include "Minimization.h"



////////////// Forward declarations to avoid large headers ////////////////
namespace pcl {

class PointXYZ;
template<class T> class PointCloud;

} // pcl


class SQParameters;


namespace ope {


/**
 * \brief Fits a parametric superquadric to a point cloud
 */
class SQFitting {
private:
	/**
	 * \brief Object to execute all minimization functions
	 */
	Minimization minProcess;


public:
	/// Constructor
	SQFitting();

	/// Destructor
	~SQFitting();
	
	
	/**
	 * \brief Initializes the Levenberg-Marquadt minimization parameters used for sq estimation
	 */
	void initializeMinimizationParameters (SQParameters& sqParams);


	/**
	 * \brief Estimate the parameters of a superquadric from the given point cloud
	 * \return the eigen vector index with the largest variance
	 */
	int estimateParameters(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams, const int& eigenVector);


	/**
	 * \brief Estimate the parameters of a superquadric from the given point cloud
	 * \return the eigen vector index with the largest variance
	 */
	int estimateParametersNew(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams, const int& eigenVector);


	/**
	 * \brief Estimate the initial SQ parameters based on eigen analysis
	 */
	void estimateInitialParameters(pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams);


	/**
	 * \brief Determines the total error (quality) of fit of superquadric parameters
	 */
	double qualityOfFit(pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& sqParams);


	/**
	 * \brief Recover the parameters of the superquadric that best fits the given point cloud
	 */
	double recoverParameters(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& prm);


	/**
	 * \brief Executes the superquadric shape fitting process on a point cloud and gets the best parameters
	 * \param <initParams> the initial superquadric parameters for the given cloud
	 * \param <bestParams> the final superquadric parameters recovered after processing
	 */
	void performShapeFitting(const pcl::PointCloud<pcl::PointXYZ>& cloud, SQParameters& initParams, SQParameters& bestParams);


};


} // ope

#endif /* SQ_FITTING_H__ */
