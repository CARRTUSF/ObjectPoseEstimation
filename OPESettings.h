/*
 * Software License Agreement (BSD License)
 *
 *  Object Pose Estimation (OPE) - www.cse.usf.edu/kkduncan
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
 *	\file	OPESettings.h
 *	\brief	Defines the settings object used throughout OPE
 *	\author	Kester Duncan
 */
#pragma once
#ifndef OPE_SETTINGS_H_
#define OPE_SETTINGS_H_


#include "OPECommon.h"

/** \namespace ope
 *	Namespace where all the Object Pose Estimation functionality resides
 */
namespace ope {

/**
 * \brief User-customizable program settings for pose estimation process
 */
struct OPESettings {
	/// Minimum depth value for target scene area. Used by the <code>PassThrough</code> filter
	float minTgtDepth;

	/// Maximum depth value for target scene area
	float maxTgtDepth;

	/// Minimum height of object hypotheses in the scene
	float minObjHeight;

	/// Maximum height of object hypotheses in the scene
	float maxObjHeight;

	/// Maximum voxel size for superquadric fitting
	float maxVoxelSize;

	/// Minimum voxel size for superquadric fitting
	float minVoxelSize;

	/// Determines whether status updates are output
	bool verbose;

	/// Determines the amount of error minimization iterations for superquadric fitting
	int minIterations;

	/// Determines whether superquadric shape tapering is allowed when estimating pose
	bool allowTapering;

	/// Object voxel size for downsampling
	float objVoxelSize;

	/// Show a debugging viewer when point cloud is captured
	bool doDebug;

	/// Default Constructor
	OPE_EXPORT OPESettings();

	/// Overloaded assignment operator
	OPE_EXPORT OPESettings& operator= (const OPESettings& other);

};

} /* ope */



#endif /* OPE_SETTINGS_H_ */