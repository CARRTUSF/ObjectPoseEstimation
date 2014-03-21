#include <cstdlib>
#include <iomanip>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/thread/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/visualization/cloud_viewer.h>

#include "TableObjectModeler.h"
#include "OPEUtils.h"

using namespace ope;


Utils::Utils() {
}


Utils::~Utils() {
}


void Utils::convertPointCloud(const pcl::PointCloud<pcl::PointXYZRGBA>& src, pcl::PointCloud<pcl::PointXYZRGB>& tgt) {
	pcl::PointXYZRGB nanPoint;
	nanPoint.getVector4fMap().setConstant (std::numeric_limits<float>::quiet_NaN());
	tgt.points.resize((int) (src.width * src.height), nanPoint);
	tgt.is_dense = src.is_dense;
	tgt.width = src.width;
	tgt.height = src.height;
	tgt.header = src.header;	

	for (size_t i = 0; i < src.points.size(); i++) {
		pcl::PointXYZRGB pt;
		pt.x = src.points[i].x;
		pt.y = src.points[i].y;
		pt.z = src.points[i].z;
		pt.r = src.points[i].r;
		pt.g = src.points[i].g;
		pt.b = src.points[i].b;

		tgt.points.push_back(pt);
	}
}


void Utils::transformPointCloud(pcl::PointCloud<pcl::PointXYZRGB>& cloud) {
	for (size_t i = 0; i < cloud.size(); ++i) {
		pcl::PointXYZRGB p = cloud.points[i];

		cloud.points[i].x = p.z;
		cloud.points[i].y = -p.x;
		cloud.points[i].z = p.y;
	}

}


void Utils::printCurrentDateTime() {
	std::ostringstream msg;
	const boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << "+-------------------------------------------+" << std::endl;
	std::cout << " " << now.date() << std::endl;
	std::cout << " " << now.time_of_day() << std::endl;
	std::cout << "+-------------------------------------------+" << std::endl << std::endl;
}


size_t Utils::getDesiredObject(pcl::PointCloud<pcl::PointXYZRGB>::ConstPtr ptCloudPtr, const std::vector<BoundingBox>& boxes) {
	srand(time(NULL));
	size_t desiredObjIdx = 0;
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Objects - In the command window, enter the object's index"));
	pcl::visualization::PointCloudColorHandlerRGBField<pcl::PointXYZRGB> rgb(ptCloudPtr);
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZRGB> (ptCloudPtr, rgb, "Cloud");
	viewer->setPointCloudRenderingProperties (pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "Cloud");

	for (size_t i = 0; i < boxes.size(); ++i) {
		std::string objId("Object_");
		std::string cubeId("Cube_");
		objId += boost::lexical_cast<std::string>(i);
		cubeId += boost::lexical_cast<std::string>(i);

		viewer->addCube(boxes[i].minX, boxes[i].maxX, boxes[i].minY, boxes[i].maxY, boxes[i].minZ, boxes[i].maxZ, 1.0, 1.0, 0.0, cubeId);
		viewer->setShapeRenderingProperties(pcl::visualization::PCL_VISUALIZER_LINE_WIDTH, 3.0, cubeId);
		pcl::PointXYZ p(boxes[i].minX, boxes[i].minY, boxes[i].minZ);
		viewer->addText3D(objId, p, 0.025);
	}
	viewer->initCameraParameters();
	viewer->resetCameraViewpoint("Cloud");

	while (!viewer->wasStopped()) {
		viewer->spinOnce(100);
		boost::this_thread::sleep (boost::posix_time::microseconds(10));
	}

	std::cout << ">> Enter the index of the object you wish to process: ";
	std::cin >> desiredObjIdx;
	while (desiredObjIdx < 0 || desiredObjIdx >= boxes.size()) {
		std::cout << "\n>> Invalid index, enter a correct object index: ";
		std::cin >> desiredObjIdx;
	}

	return desiredObjIdx;

}