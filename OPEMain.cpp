/**
 * \file Main.cpp
 * \brief Main entry point for Object Pose Estimation
 * \author Kester Duncan
 *
 * This file is the main entry point for using the object pose estimator using superquadrics by 
 * Kester Duncan
 */
#if 1
#include <iostream>
#include "SQTypes.h"
#include "ObjectPoseEstimator.h"


int main (int argc, char *argv[]) {
	ope::OPESettings settings;
	settings.minTgtDepth = 0.4f;
	settings.maxObjHeight = 2.5f;

	ope::ObjectPoseEstimator estimator(settings);
	ope::SQParameters sqParams = estimator.run();
		
	std::cin.get();
	std::cin.get();

	return 0;
}


#endif