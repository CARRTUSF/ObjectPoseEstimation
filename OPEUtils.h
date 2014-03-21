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
 *	\file Utils.h
 *	\author Kester Duncan
 */
#pragma once
#ifndef OPE_UTILS_H__
#define OPE_UTILS_H__


////////////// Forward declarations to avoid large headers ////////////////
namespace pcl {
	class PointXYZ;
	class PointXYZRGB;
	class PointXYZRGBA;
	template<class T> class PointCloud;

} // pcl



namespace ope {

// Forward declaration to avoid including header
class BoundingBox;


/**
 * \brief Performs frequently used utility functions
 * \author Kester Duncan
 */
class Utils
{
public:
	Utils();
	~Utils();


	/**
	 * \brief Convert between pcl point cloud types
	 * \details Converts from pcl::PointXYZRGBA to pcl::PointXYZRGB
	 * \param <src> - input <code>pcl::PointXYZRGBA</code> point cloud
	 * \param <tgt> - output <code>pcl::PointXYZRGB</code> point cloud
	 */
	static void convertPointCloud(const pcl::PointCloud<pcl::PointXYZRGBA>& src, pcl::PointCloud<pcl::PointXYZRGB>& tgt); 


	/**
	 * \brief Transform an XYZRGB point cloud from the Kinect optical frame to world coordinate frame
	 * \param <cloud> - the cloud to be transformed
	 */
	static void transformPointCloud(pcl::PointCloud<pcl::PointXYZRGB>& cloud);


	/**
	 * \brief Prints the current local time to the output stream
	 */
	static void printCurrentDateTime();


	/**
	 * \brief Displays a pointcloud with highlighted objects in order to determine the user's object of choice
	 */
	static size_t getDesiredObject(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr ptCloudPtr, const std::vector<BoundingBox>& boxes);

};


} /* ope */


#endif /* OPE_UTILS_H__ */

