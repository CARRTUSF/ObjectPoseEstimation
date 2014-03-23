#include <iostream>

#include <boost/progress.hpp>

#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types.h>
#include <pcl/surface/mls.h>
#include <pcl/surface/impl/mls.hpp>

#include "OPESettings.h"
#include "SQTypes.h"
#include "OPEUtils.h"
#include "Minimization.h"
#include "SQFitting.h"
#include "PointCloudCapture.h"
#include "TableObjectModeler.h"
#include "ObjectPoseEstimator.h"


namespace ope {

ObjectPoseEstimator::ObjectPoseEstimator() : settings(OPESettings()) {
	
}


ObjectPoseEstimator::ObjectPoseEstimator(const OPESettings& settings) {
	this->settings = settings;

}


SQParameters ObjectPoseEstimator::calculateObjectPose(pcl::PointCloud<pcl::PointXYZRGB>& selectedObjectPtCloud) {
	// Main object instance for superquadric fitting
	SQFitting fittingProcess;

	// The initial XYZ point cloud that is processes for pose estimation
	pcl::PointCloud<pcl::PointXYZ> cloud;

	// The downsampled XYZ point cloud that is used in the multiscale voxelization scheme
	pcl::PointCloud<pcl::PointXYZ> downsampledCloud;

	// A pointer to the initial XYZ point cloud used during voxelization
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloudPtr (new pcl::PointCloud<pcl::PointXYZ>());

	// Initial superquadric parameters based on the cloud properties
	SQParameters initParams;	

	// The final (best) superquadric parameters that are found
	SQParameters bestParams;

	// Estimate of the superquadric parameters used for bootstrapping the multiscale voxelization scheme
	SQParameters sqEstimate;

	// A record of the initial number of points before repeated downsampling
	int numberOfPoints = 0;

	/*
	 * The point cloud must be transformed from the Kinect optical frame to the
	 * world coordinate frame. Additional transformations may have to be done depending
	 * on the application, but this is left up to the user.
	 */
	Utils::transformPointCloud(selectedObjectPtCloud);

	/*
	 * The captured XYZRGB point cloud must be converted to an XYZ cloud
	 * in order to perform pose estimation
	 */
	for (size_t i = 0; i < selectedObjectPtCloud.size(); ++i) {
		pcl::PointXYZ p;
		p.x = selectedObjectPtCloud[i].x;
		p.y = selectedObjectPtCloud[i].y;
		p.z = selectedObjectPtCloud[i].z;

		cloud.points.push_back(p);
		numberOfPoints++;
	}	
	cloud.width = numberOfPoints;
	cloud.height = 1;
	*cloudPtr = cloud;

	/*
	 * The voxel downsampling parameters must be initialized
	 */
	const int NUM_SCALES = 4; // 4 is our set number
	const float MAX_VOXEL = settings.maxVoxelSize;
	const float MIN_VOXEL = settings.minVoxelSize;
	const float SCALE_DIFF = (MAX_VOXEL - MIN_VOXEL) / (float) NUM_SCALES;
	float gridSizes[NUM_SCALES];
	double totalElapsedTime = 0;

	for (int s = 0; s < NUM_SCALES; s++) {
		gridSizes[s] = MAX_VOXEL - (s * SCALE_DIFF);
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr smoothCloudPtr (new pcl::PointCloud<pcl::PointXYZ>);
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
	pcl::MovingLeastSquares<pcl::PointXYZ, pcl::PointXYZ> mls;
	mls.setComputeNormals (false);

	// Set parameters
	mls.setInputCloud (cloudPtr);
	mls.setPolynomialFit (true);
	mls.setSearchMethod (tree);
	mls.setSearchRadius (0.03);

	// Reconstruct
	mls.process (*smoothCloudPtr);

	/* 
	 * Calculate first error differences
	 */
	/// Variables pertaining to multi-scale voxelization algorithm
	double errorThreshold = 2.0;
	double errorDiff = 1000.0;
	double errorValue = 0;
	double prevErrorValue = 0;

	fittingProcess.estimateInitialParameters(*smoothCloudPtr, sqEstimate);
	errorValue = fittingProcess.qualityOfFit(*smoothCloudPtr, sqEstimate);
	errorDiff = abs(prevErrorValue - errorValue);
	prevErrorValue = errorValue;

	if (settings.verbose == true) {
		std::cout << ">> Performing object pose estimation using Superquadrics\n";
	}

	/*
	 * Multi-scale voxelization for pose estimation
	 */
	for (int j = 0; j < NUM_SCALES && errorDiff >= errorThreshold; j++) {
		if (j != NUM_SCALES - 1) {
			pcl::VoxelGrid<pcl::PointXYZ> grid;
			grid.setInputCloud (smoothCloudPtr);
			grid.setLeafSize (gridSizes[j], gridSizes[j], gridSizes[j]);
			grid.filter (downsampledCloud);

		} else {
			downsampledCloud = *smoothCloudPtr;
		}

		/*
		 * Estimate the dimensions of cloud (in order to get a more accurate fit)
		 */
		fittingProcess.estimateInitialParameters(downsampledCloud, sqEstimate);
		
		/*
		 * Initialize the superquadric parameters based on the cloud dimensions
		 */
		initParams.min.iterations = settings.minIterations;		
		initParams.min.a1.type = BOUNDED;
		initParams.min.a1.lowerBound = 0.020f;
		initParams.min.a1.upperBound = sqEstimate.a1 + 0.015f;
		initParams.a1 = 0.05f;
	
		initParams.min.a2.type = BOUNDED;
		initParams.min.a2.lowerBound = 0.020f;
		initParams.min.a2.upperBound = sqEstimate.a2 + 0.015f;
		initParams.a2 = 0.05f;

		initParams.min.a3.type = BOUNDED;
		initParams.min.a3.lowerBound = 0.020f;
		initParams.min.a3.upperBound = sqEstimate.a3 + 0.015f;
		initParams.a3 = 0.05f;

		initParams.min.e1.type = BOUNDED;
		initParams.min.e1.lowerBound = 0.1f;
		initParams.min.e1.upperBound = 1.0f;
		initParams.e1 = 1.0f;

		initParams.min.e2.type = BOUNDED;
		initParams.min.e2.lowerBound = 0.1f;
		initParams.min.e2.upperBound = 1.0f;
		initParams.e2 = 1.0f;
	
		initParams.min.phi.type = UNLIMITED;
		initParams.min.theta.type = UNLIMITED;
		initParams.min.psi.type = UNLIMITED;
		initParams.min.phi.value = 1.0f;
		initParams.min.theta.value = 1.0f;
		initParams.min.psi.value = 1.0f;
	
		initParams.min.px.type = UNLIMITED;
		initParams.min.py.type = UNLIMITED;
		initParams.min.pz.type = UNLIMITED;
		
		if (settings.allowTapering) {
			initParams.min.kx.type = BOUNDED;
			initParams.min.kx.lowerBound = -0.25f;
			initParams.min.kx.upperBound = 0.25f;
			initParams.min.kx.value = 0.0f;
			initParams.kx = 1.0f;
		
			initParams.min.ky.type = BOUNDED;
			initParams.min.ky.lowerBound = -0.25f;
			initParams.min.ky.upperBound = 0.25;
			initParams.min.ky.value = 0.0f;
			initParams.ky = 1.0f;

		} else {
			initParams.min.kx.type = UNCHANGED;
			initParams.min.kx.value = 0.0f;	
			initParams.min.ky.type = UNCHANGED;
			initParams.min.ky.value = 0.0f;

		}		

		std::stringstream os;
		double elapsedTime = 0.0;            
		{
			boost::progress_timer timer(os);
			fittingProcess.performShapeFitting(downsampledCloud, initParams, bestParams);
		}

		if (!(os >> elapsedTime)) {
			elapsedTime = atof(os.str().c_str());
		}

		totalElapsedTime += elapsedTime;
			
		errorValue = fittingProcess.qualityOfFit(cloud, bestParams);
		errorDiff = abs(prevErrorValue - errorValue);
		prevErrorValue = errorValue;

		/*
		 * The parameters found at this level are used to initialize the successive one
		 */
		bestParams.copyTo(initParams);
		downsampledCloud.clear();
	}

	initParams.copyTo(bestParams);	

	if (settings.verbose == true) {
		std::cout <<">> Superquadric fitting processing time: " << totalElapsedTime << " seconds" << std::endl;
	}
	
	return bestParams;
}


SQParameters ObjectPoseEstimator::run() {
	SQParameters sqParams;
	PointCloudCapture cap;
	pcl::PointCloud<pcl::PointXYZRGB>::Ptr cloudPtr (new pcl::PointCloud<pcl::PointXYZRGB>);
	bool detectedObjects;
	size_t desiredObjIdx = 0;

	/*
	 * Capture point cloud
	 */
	cap.run(*cloudPtr);
	

	/*
	 * Generate object models based on extracted clusters
	 */
	ObjectModelGenerator<pcl::PointXYZRGB> generator(cloudPtr);
	detectedObjects = generator.generateObjectModels(settings);


	/*
	 * Choose the object whose pose (position and orientation) you wish to calculate
	 */
	if (detectedObjects) {
		desiredObjIdx = Utils::getDesiredObject(cloudPtr, generator.getBoundingBoxes());
		
		/*
		 * Calculate the object pose using Superquadrics
		 */
		sqParams = ObjectPoseEstimator::calculateObjectPose(generator.objects[desiredObjIdx].objectCloud);
		std::cout << sqParams;

	} else {
		if (settings.verbose) {
			std::cout << ">> NO OBJECTS WERE DETECTED!" << std::endl;
		}
	}

	return sqParams;
}


} // ope