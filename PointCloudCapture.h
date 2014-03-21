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
 *	\file	PointCloudCapture.h
 *	\author	Kester Duncan
 */
#pragma once
#ifndef POINTCLOUDCAPTURE_H__
#define POINTCLOUDCAPTURE_H__

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

namespace ope {

/**
 * \brief Captures XYZRGB point clouds from the Kinect
 * \note The coordinate system for the captured point clouds is as follows:
 * x-axis -> right, y-axis -> down, z-axis -> points into scene.
 *                          ________+x
 *                         |\
 *                         | \ 
 *                         |  \
 *                         |   \
 *                         |    \
 *                         +y    +z
 *
 * Also, a PassThrough filter is applied to the captured point cloud so that the 
 * target area is within a specified range which is adjustable
 */
class PointCloudCapture
{
private:
	pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr ptCloudPtr;
	pcl::PointCloud<pcl::PointXYZRGBA>::Ptr filteredPtCloudPtr;
	
public:
	PointCloudCapture();
	~PointCloudCapture();
	
	/**
	 * \brief Captures an XYZRGBA point cloud from the Kinect and stores an XYZRGB cloud
	 * \param <ptCloud> - holds the captured PointXYZRGB point cloud
	 */
	void run(pcl::PointCloud<pcl::PointXYZRGB>& ptCloud, const OPESettings& settings = OPESettings());


	/// \brief Callback function that is used to read XYZRGB point clouds from the Kinect
	void cloudCallback(const pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr &cloud);
};


} /* end of namespace ope */


#endif /* POINTCLOUDCAPTURE_H__ */

