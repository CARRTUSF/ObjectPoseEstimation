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
 *	\file	SQTypes.h
 *	\brief	Defines types used for Superquadric shape fitting
 *	\author	Kester Duncan
 */
#pragma once
#ifndef __SHAPEFITTINGTYPES_H__
#define __SHAPEFITTINGTYPES_H__

#include <iostream>
#include "OPESettings.h"

/** \namespace ope
 *	Namespace where all the Object Pose Estimation functionality resides
 */
namespace ope {

/// Defines the numerical property of a superquadric parameter
typedef enum ParameterType {
	UNCHANGED = 1,
	BOUNDED,
	UNLIMITED,
	NOT_USED,
	OFFSET

} ParameterType;


/// Defines the upper and lower limits for superquadric parameters
typedef struct ParameterLimits {
	ParameterType type;			
	double value;
	double lowerBound;
	double upperBound;
	ParameterLimits() : value(0.f), lowerBound(0.f), upperBound(0.f) {}

} ParameterLimits;


/// Defines the minimization parameters that would be used for minimization during SQ fitting
typedef struct MinimizationParameters {
	int iterations;
	ParameterLimits a1;
	ParameterLimits a2;
	ParameterLimits a3;
	
	ParameterLimits e1;
	ParameterLimits e2;

	ParameterLimits phi;
	ParameterLimits theta;
	ParameterLimits psi;

	ParameterLimits px;
	ParameterLimits py;
	ParameterLimits pz;

	ParameterLimits kx;
	ParameterLimits ky;

	MinimizationParameters() : iterations(20) {}

} MinimizationParameters;


/**
 * \brief Defines the Superquadric parameters used for object pose estimation
 *
 */
class SQParameters{
public:
	/// The shape dimension for the x-axis
	double a1;
	/// The shape dimension for the y-axis
	double a2;
	/// The shape dimension for the z-axis
	double a3;

	/// The north-south superquadric shape parameter
	double e1;
	/// The east-west superquadric shape parameter
	double e2;

	/// The x-axis location of the centroid of this superquadric
	double px;
	/// The y-axis location of the centroid of this superquadric
	double py;
	/// The z-axis location of the centroid of this superquadric
	double pz;

	/// Euler rotation angle along the x-axis
	double phi;
	/// Euler rotation angle along the y-axis
	double theta;
	/// Euler rotation angle along the x-axis
	double psi;

	/// Tapering parameter along the x-axis
	double kx;
	/// Tapering parameter along the y-axis
	double ky;

	/// Index of the principal axis in the calculated rotation matrix
	int principalAxis;
	int majorAxis;
	int minorAxis;

	/// Minimization properties
	MinimizationParameters min;


	/// Default Constructor
	SQParameters() {		
		a1 = 0.0; a2 = 0.0; a3 = 0.0;
		e1 = 1.0; e2 = 1.0; 
		px = 0.0; py = 0.0; pz = 0.0;
		phi = 0.0; theta = 0.0; psi = 0.0;
		kx = 0.0; ky = 0.0;
		principalAxis = 0;
		majorAxis = 0;
		minorAxis = 0;

		min.iterations = 20;
		min.a1.type = BOUNDED;
		min.a1.lowerBound = 0.020;
		min.a1.upperBound = 0.30;
		a1 = 0.05;
	
		min.a2.type = BOUNDED;
		min.a2.lowerBound = 0.020;
		min.a2.upperBound = 0.30;
		a2 = 0.05;

		min.a3.type = BOUNDED;
		min.a3.lowerBound = 0.020;
		min.a3.upperBound = 0.30;
		a3 = 0.05;

		min.e1.type = BOUNDED;
		min.e1.lowerBound = 0.1;
		min.e1.upperBound = 2.0;
		e1 = 1;

		min.e2.type = BOUNDED;
		min.e2.lowerBound = 0.1f;
		min.e2.upperBound = 2.0f;
		e2 = 1;
	
		min.phi.type = UNLIMITED;
		min.theta.type = UNLIMITED;
		min.psi.type = UNLIMITED;
		min.phi.value = 1.0;
		min.theta.value = 1.0;
		min.psi.value = 1.0;
	
		min.px.type = UNLIMITED;
		min.py.type = UNLIMITED;
		min.pz.type = UNLIMITED;		
		
		min.kx.type = BOUNDED;
		min.kx.lowerBound = -1.0;
		min.kx.upperBound = 1.0;
		min.kx.value = 0.0;
		kx = 1.0;

		min.ky.type = BOUNDED;
		min.ky.lowerBound = -1.0;
		min.ky.upperBound = 1.0;
		min.ky.value = 0.0;
		ky = 1.0;

		
	
	}


	/// Deep copy of this class' properties
	void copyTo(SQParameters& sqParams) {
		sqParams.a1 = this->a1;
		sqParams.a2 = this->a2;
		sqParams.a3 = this->a3;
		sqParams.e1 = this->e1;
		sqParams.e2 = this->e2;
		sqParams.psi = this->psi;
		sqParams.theta = this->theta;
		sqParams.phi = this->phi;
		sqParams.px = this->px;
		sqParams.py = this->py;
		sqParams.pz = this->pz;
		sqParams.kx = this->kx;
		sqParams.ky = this->ky;
		sqParams.principalAxis = this->principalAxis;
		sqParams.majorAxis = this->majorAxis;
		sqParams.minorAxis = this->minorAxis;

		sqParams.min.a1.lowerBound = this->min.a1.lowerBound;
		sqParams.min.a1.upperBound = this->min.a1.upperBound;
		sqParams.min.a1.type = this->min.a1.type;
		sqParams.min.a1.value = this->min.a1.value;

		sqParams.min.a2.lowerBound = this->min.a2.lowerBound;
		sqParams.min.a2.upperBound = this->min.a2.upperBound;
		sqParams.min.a2.type = this->min.a2.type;
		sqParams.min.a2.value = this->min.a2.value;

		sqParams.min.a3.lowerBound = this->min.a3.lowerBound;
		sqParams.min.a3.upperBound = this->min.a3.upperBound;
		sqParams.min.a3.type = this->min.a3.type;
		sqParams.min.a3.value = this->min.a3.value;

		sqParams.min.e1.lowerBound = this->min.e1.lowerBound;
		sqParams.min.e1.upperBound = this->min.e1.upperBound;
		sqParams.min.e1.type = this->min.e1.type;
		sqParams.min.e1.value = this->min.e1.value;

		sqParams.min.e2.lowerBound = this->min.e2.lowerBound;
		sqParams.min.e2.upperBound = this->min.e2.upperBound;
		sqParams.min.e2.type = this->min.e2.type;
		sqParams.min.e2.value = this->min.e2.value;

		sqParams.min.psi.lowerBound = this->min.psi.lowerBound;
		sqParams.min.psi.upperBound = this->min.psi.upperBound;
		sqParams.min.psi.type = this->min.psi.type;
		sqParams.min.psi.value = this->min.psi.value;

		sqParams.min.theta.lowerBound = this->min.theta.lowerBound;
		sqParams.min.theta.upperBound = this->min.theta.upperBound;
		sqParams.min.theta.type = this->min.theta.type;
		sqParams.min.theta.value = this->min.theta.value;

		sqParams.min.phi.lowerBound = this->min.phi.lowerBound;
		sqParams.min.phi.upperBound = this->min.phi.upperBound;
		sqParams.min.phi.type = this->min.phi.type;
		sqParams.min.phi.value = this->min.phi.value;

		sqParams.min.px.lowerBound = this->min.px.lowerBound;
		sqParams.min.px.upperBound = this->min.px.upperBound;
		sqParams.min.px.type = this->min.px.type;
		sqParams.min.px.value = this->min.px.value;

		sqParams.min.py.lowerBound = this->min.py.lowerBound;
		sqParams.min.py.upperBound = this->min.py.upperBound;
		sqParams.min.py.type = this->min.py.type;
		sqParams.min.py.value = this->min.py.value;

		sqParams.min.pz.lowerBound = this->min.pz.lowerBound;
		sqParams.min.pz.upperBound = this->min.pz.upperBound;
		sqParams.min.pz.type = this->min.pz.type;
		sqParams.min.pz.value = this->min.pz.value;

		sqParams.min.kx.lowerBound = this->min.kx.lowerBound;
		sqParams.min.kx.upperBound = this->min.kx.upperBound;
		sqParams.min.kx.type = this->min.kx.type;
		sqParams.min.kx.value = this->min.kx.value;

		sqParams.min.ky.lowerBound = this->min.ky.lowerBound;
		sqParams.min.ky.upperBound = this->min.ky.upperBound;
		sqParams.min.ky.type = this->min.ky.type;
		sqParams.min.ky.value = this->min.ky.value;
		sqParams.min.iterations = this->min.iterations;
	}


	/**
	 * \brief Prints the superquadric parameters to an output stream
	 */
	friend std::ostream& operator<<(std::ostream &os, const SQParameters& sq) {
		os << ">>\t Superquadric Parameters:\n";
		os << "\t   [location]    (px, py, pz) --> (" << sq.px << ", " << sq.py << ", " << sq.pz << ")\n";
		os << "\t   [dimensions]  (a1, a2, a3) --> (" << sq.a1 << ", " << sq.a2 << ", " << sq.a3 << ")\n";
		os << "\t   [rotation]    (phi, theta, psi) --> (" << sq.phi << ", " << sq.theta << ", " << sq.psi << ")\n";
		os << "\t   [shape]       (e1, e2) --> (" << sq.e1 << ", " << sq.e2 << ")\n";
		os << "\t   [tapering]    (kx, ky) --> (" << sq.kx << ", " << sq.ky << ")\n";		
		os << std::endl;

		return os;
	}

};


} // ope

#endif /* __SHAPEFITTINGTYPES_H__ */