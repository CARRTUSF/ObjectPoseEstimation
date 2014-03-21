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
 *	\file	TableObjectModeler.h
 *	\brief	Uses the output of the TableObjectDetector to create generic object models
 *	\author	Kester Duncan
 */
#pragma once
#ifndef __TABLE_OBJECT_MODELER_H__
#define __TABLE_OBJECT_MODELER_H__

#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/filters/statistical_outlier_removal.h>
#include <boost/lexical_cast.hpp>
#include "TableObjectDetector.h"


namespace ope {

/**
 * \brief Defines a 3D bounding box around a point cloud
 */
class BoundingBox {
public:
	BoundingBox() : minX(999), maxX(-999), minY(999), maxY(-999), minZ(999), maxZ(-999), width(0), height(0), depth(0) {}
	BoundingBox(float _minX, float _maxX, float _minY, float _maxY, float _minZ, float _maxZ, float _w, float _h, float _d)
		: minX(_minX), maxX(_maxX), minY(_minY), maxY(_maxY), minZ(_minZ), maxZ(_maxZ), width(_w), height(_h), depth(_d) {}
	
	
	inline bool isEmpty() const {
        return (width < 0 || height < 0 || depth < 0);
    }


	inline pcl::PointXYZ centroid() const {
        return pcl::PointXYZ((minX + maxX) / 2.0f, (minY + maxY) / 2.0f, (minZ + maxZ) / 2.0f);
    }


    inline bool isPointInside(const pcl::PointXYZ& p) const {
        return ((p.x >= minX) && (p.x <= maxX) && (p.y >= minY) && (p.y <= maxY) && (p.z >= minZ) && (p.z <= maxZ));
    }

	inline bool isPointInside2D(const int& x, const int& y) {
		bool isInside = false;
		int pixelOffset = 5;
		
		if (x >= (xPos() - pixelOffset) && x <= (xPos() + pixelOffset) && 
			y >= (yPos() - pixelOffset) && y <= (yPos() + pixelOffset)) {
			isInside = true;
		}
		return isInside;
	}


	inline float volume() const {
		return (width * height * depth);
	}


	inline void update(float x, float y, float z) {
		if (x < minX) minX = x;
		if (y < minY) minY = y;
		if (z < minZ) minZ = z;

		if (x > maxX) maxX = x;
		if (y > maxY) maxY = y;
		if (z > maxZ) maxZ = z;

		width = maxX - minX;
		height = maxY - minY;
		depth = maxZ - minZ;
	}


	inline int xPos() const {
		// 2D screen position
		float xPos = (525.0f * minX) / minZ;
		xPos += 320;
		
		return ((int) xPos);
	}


	inline int yPos() const {
		// 2D screen position
		float yPos = (525.0f * minY) / minZ;
		yPos += 240;		

		return ((int) yPos);
	}


	friend std::ostream& operator<<(std::ostream &os, const BoundingBox& b) {
		os << "\n\tBounding Box Properties:\n";
		os << "\t\t (min X, max X) --> (" << b.minX << ", " << b.maxX << ")\n";
		os << "\t\t (min Y, max Y) --> (" << b.minY << ", " << b.maxY << ")\n";
		os << "\t\t (min Z, max Z) --> (" << b.minZ << ", " << b.maxZ << ")\n";
		os << "\t\t (width, height, depth) --> (" << b.width << ", " << b.height << ", " << b.depth << ")\n";
		os << "\t\t (x position, y position) --> (" << b.xPos() << ", " << b.yPos() << ")\n";		
		os << std::endl;

		return os;
	}

	
public:
    float minX, maxX;
	float minY, maxY;
	float minZ, maxZ;
	float width;
	float height;
	float depth;

}; // BoundingBox




/**
 * \brief Models object hypotheses found on a table
 */
template <class PointType>
class TableObjectModel {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	pcl::PointCloud<PointType> objectCloud;
	BoundingBox box;
	size_t objectId;	///< A cluster object's unique id according to the order it was detected
	std::string name;	///< Cluster's identifier

public:
	TableObjectModel() : objectId(999), name("Unknown") {}
	~TableObjectModel() {
		if (!objectCloud.empty()) {
			objectCloud.clear();			
		}	
	}
	

	/**
	 * \brief Create the object model
	 */
	void create(const TableObjectDetector< PointType >& detector, const size_t& clusterId) {
		if (clusterId >= detector.objectClusters().size()) {
			std::cerr << "TableObjectModeler: Invalid cluster!" << std::endl;
		}

		pcl::PointCloud< PointType >::Ptr tempCloudPtr (new pcl::PointCloud< PointType>);
				
		objectId = clusterId;
		name = "Object_";
		name += boost::lexical_cast<std::string>(objectId);
		
		size_t numPoints = detector.objectClusters()[clusterId].size();
		tempCloudPtr->reserve(numPoints);

		for (size_t i = 0; i < numPoints; ++i) {
			PointType p = detector.objectClusters()[clusterId][i];						
			tempCloudPtr->push_back(p);
		}

		pcl::StatisticalOutlierRemoval< PointType > sor (true);
		sor.setInputCloud (tempCloudPtr);
		sor.setMeanK(8);
		sor.setStddevMulThresh(1.0);
		sor.filter(objectCloud);

		for (size_t i = 0; i < objectCloud.size(); ++i) {
			PointType p = objectCloud.points[i];
			box.update(p.x, p.y, p.z);			
		}		
	}



}; // TableObjectModel



/**
 * \brief Generate point cloud models for objects detected on a table
 */
template <class PointType>
class ObjectModelGenerator {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	typedef typename pcl::PointCloud<PointType> PointCloud;
	typedef typename PointCloud::Ptr PointCloudPtr;
	typedef typename PointCloud::ConstPtr PointCloudConstPtr;

private:
	PointCloudConstPtr srcCloudPtr;
	
public:
	/**
	 * \brief List of object models
	 *
	 * \note This declaration is rather messy but it is required if we desire to use ANYTHING Eigen. All
	 * Eigen objects are 16-byte aligned. Therefore we must use Eigen's memory allocator when using standard
	 * containers.
	 */
	std::vector < TableObjectModel< PointType >, Eigen::aligned_allocator< TableObjectModel< PointType > > > objects;


	ObjectModelGenerator() : detectedObjects(false) { }
	ObjectModelGenerator(PointCloudConstPtr ptCloudPtr) {
		setSourceCloud(ptCloudPtr);		
	}

	~ObjectModelGenerator() {
		if(!objects.empty()) {
			objects.clear();
		}
	}
	

	void setSourceCloud(PointCloudConstPtr ptCloudPtr) {
		srcCloudPtr = ptCloudPtr;
		
	}


	bool generateObjectModels(const OPESettings& settings = OPESettings()) {
		TableObjectDetector< PointType > detector;
		bool detectedObjects = detector.detect(srcCloudPtr);

		if (settings.verbose) {
			std::cout << ">> Generating object models" << std::endl;

			if (detectedObjects) {
				std::cout << ">>\t " << detector.objectClusters().size() << " object(s) detected" << endl;
			} else {
				std::cout << ">>\t " << "No objects were detected" << std::endl;
			}
		}

		for (size_t i = 0; i < detector.objectClusters().size(); ++i) {
			TableObjectModel< PointType > objectModel;
			objectModel.create(detector, i);
			objects.push_back(objectModel);
		}

		if (detector.objectClusters().size() > 0) {
			detectedObjects = true;
		} else {
			detectedObjects = false;
		}

		return detectedObjects;
	}


	std::vector<BoundingBox> getBoundingBoxes() {
		std::vector<BoundingBox> objectBoxes;

		for (size_t i = 0; i < objects.size(); ++i) {
			objectBoxes.push_back(objects[i].box);
		}

		return objectBoxes;
	}
	
}; // ObjectModelGenerator



} /* ope */


#endif /* __TABLE_OBJECT_MODELER_H__ */