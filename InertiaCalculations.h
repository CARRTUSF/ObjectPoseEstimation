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
 *	\file	InertiaCalculations.h
 *	\brief	Functions to calculate the moments of inertia of superquadrics. 
 *			This is used to calculate the principal axis of a superquadric, which 
 *			coincides with the axis of least inertia
 *	\author	Kester Duncan
 */
#pragma once
#ifndef _INERTIA_CALCULATIONS_H_
#define _INERTIA_CALCULATIONS_H_


namespace ope {


// Forward declaration
class SQParameters;


struct InertiaCalculations {
	/**
	 * \brief Beta function (or Euler integral of the first kind)
	 * For more information, visit: http://en.wikipedia.org/wiki/Beta_function
	 */
	double betaFunction (const double& x, const double& y);
	

	/**
	 * \brief Calculates the moment of inertia along the x axis
	 */
	double xMomentOfInertia(SQParameters& sq);


	/**
	 * \brief Calculates the moment of inertia along the y axis
	 */
	double yMomentOfInertia(SQParameters& sq);


	/**
	 * \brief Calculates the moment of inertia along the z axis
	 */
	double zMomentOfInertia(SQParameters& sq);


	/**
	 * \brief Calculates the principal axis of a superquadric
	 */
	void calculateNewPrincipalAxis (SQParameters& sq);


	/**
	 * \brief Calculates the major axis of a superquadric
	 */
	void calculateMajorAxis (SQParameters& sq);


	/**
	 * \brief Calculates the minor axis of a superquadric
	 */
	void calculateMinorAxis (SQParameters& sq);


};


} /* ope namespace */

#endif /*_INERTIA_CALCULATIONS_H_ */