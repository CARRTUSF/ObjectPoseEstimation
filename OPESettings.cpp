#include "OPESettings.h"

namespace ope {

OPESettings::OPESettings() : minTgtDepth(0.2f), maxTgtDepth(1.8f), minObjHeight(0.01f), maxObjHeight(0.50f),
	objVoxelSize(0.003f), maxVoxelSize(0.03f), minVoxelSize(0.003f), 
	verbose(false), allowTapering(false), doDebug(false), minIterations(20) {

}


OPESettings& OPESettings::operator= (const OPESettings& other) {
	this->allowTapering = other.allowTapering;
	this->doDebug = other.doDebug;
	this->verbose = other.verbose;
	this->minObjHeight = other.minObjHeight;
	this->maxObjHeight = other.maxObjHeight;
	this->minTgtDepth = other.minTgtDepth;
	this->maxTgtDepth = other.maxTgtDepth;
	this->minVoxelSize = other.minVoxelSize;
	this->maxVoxelSize = other.maxVoxelSize;
	this->minIterations = other.minIterations;	
	this->objVoxelSize = other.objVoxelSize;

	return *this;

}


}