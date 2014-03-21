/**
 * \file TableObjectDetector.hpp
 * \author Nicolas Burrus, altered by Kester Duncan
 * \note Adapted from ntk's TableObjectDetector
 *
 * \details This implementation file is used to separate the declared templates from
 * the implementation but at the same time, make the implementation of the templates
 * visible to the compiler. This is so because if the implementation is done in a
 * *.cpp file, the compiler cannot explicitly instantiate a new class given the
 * template argument. The implementation must be found with the declaration, which
 * is a stupid but necessary C++ limitation.
 */
#include <boost/make_shared.hpp>
#include "TableObjectDetector.h"


namespace ope {

template <class PointType>
TableObjectDetector<PointType> :: TableObjectDetector() : m_max_dist_to_plane(0.03), settings(OPESettings()) {
	
	// ---[ Create all PCL objects and set their parameters
	setObjectVoxelSize(settings.objVoxelSize);
	setBackgroundVoxelSize();
	setDepthLimits(settings.minTgtDepth, settings.maxTgtDepth);
	setObjectHeightLimits(settings.minObjHeight, settings.maxObjHeight);

	// Normal estimation parameters
	k_ = 15; // 50 k-neighbors by default

	// Table model fitting parameters
	sac_distance_threshold_ = 0.01; // 1cm

	// Don't know ?
	normal_distance_weight_ = 0.1;

	// Clustering parameters
	object_cluster_tolerance_ = 0.05;	// 5cm between two objects
	object_cluster_min_size_  = 100;	// 100 points per object cluster
}


template <class PointType>
void TableObjectDetector<PointType> :: initialize() {
	grid_.setLeafSize (downsample_leaf_, downsample_leaf_, downsample_leaf_);
	grid_objects_.setLeafSize (downsample_leaf_objects_, downsample_leaf_objects_, downsample_leaf_objects_);
	grid_.setFilterFieldName ("z");
	pass_.setFilterFieldName ("z");

	grid_.setFilterLimits (min_z_bounds_, max_z_bounds_);
	pass_.setFilterLimits (min_z_bounds_, max_z_bounds_);
	grid_.setDownsampleAllData (false);
	grid_objects_.setDownsampleAllData (false);

	// Only works if PCL version > 1.2.0
	normals_tree_ = boost::make_shared<pcl::search::KdTree<Point> > ();
	clusters_tree_ = boost::make_shared<pcl::search::KdTree<Point> > ();
	clusters_tree_->setEpsilon (1);	

	n3d_.setKSearch (k_);
	n3d_.setSearchMethod (normals_tree_);

	normal_distance_weight_ = 0.1;

	seg_.setNormalDistanceWeight (normal_distance_weight_);
	seg_.setOptimizeCoefficients (true);
	seg_.setModelType (pcl::SACMODEL_NORMAL_PLANE);
	seg_.setMethodType (pcl::SAC_RANSAC);
	seg_.setProbability (0.99);
	seg_.setDistanceThreshold (sac_distance_threshold_);
	seg_.setMaxIterations (10000);

	proj_.setModelType (pcl::SACMODEL_NORMAL_PLANE);

	prism_.setHeightLimits (object_min_height_, object_max_height_);

	cluster_.setClusterTolerance (object_cluster_tolerance_);
	cluster_.setMinClusterSize (object_cluster_min_size_);
	cluster_.setSearchMethod (clusters_tree_);
}


template <class PointType>
bool TableObjectDetector<PointType> :: detect(PointCloudConstPtr cloud) {
	m_object_clusters.clear();
	initialize();

	// ---[ Convert the dataset
	cloud_ = cloud; 

	// ---[ Create the voxel grid
	pcl::PointCloud<Point> cloud_filtered;
	pass_.setInputCloud (cloud_);
	pass_.filter (cloud_filtered);
	cloud_filtered_.reset (new pcl::PointCloud<Point> (cloud_filtered));
	
	pcl::PointCloud<Point> cloud_downsampled;
	grid_.setInputCloud (cloud_filtered_);
	grid_.filter (cloud_downsampled);
	cloud_downsampled_.reset (new pcl::PointCloud<Point> (cloud_downsampled));

	if ((int) cloud_filtered_->points.size () < k_) {
		// This means that filtering returned very few points
		return false;
	}

	// ---[ Estimate the point normals
	pcl::PointCloud<pcl::Normal> cloud_normals;
	n3d_.setInputCloud (cloud_downsampled_);
	n3d_.compute (cloud_normals);
	cloud_normals_.reset (new pcl::PointCloud<pcl::Normal> (cloud_normals));
	
	// ---[ Perform segmentation
	pcl::PointIndices table_inliers; pcl::ModelCoefficients table_coefficients;
	seg_.setInputCloud (cloud_downsampled_);
	seg_.setInputNormals (cloud_normals_);
	seg_.segment (table_inliers, table_coefficients);
	table_inliers_.reset (new pcl::PointIndices (table_inliers));
	table_coefficients_.reset (new pcl::ModelCoefficients (table_coefficients));
	
	if (table_inliers_->indices.size () == 0)
		return false;

	m_plane = Plane (table_coefficients.values[0],
		table_coefficients.values[1],
		table_coefficients.values[2],
		table_coefficients.values[3]);

	// ---[ Extract the table
	pcl::PointCloud<Point> table_projected;
	proj_.setInputCloud (cloud_downsampled_);
	proj_.setIndices (table_inliers_);
	proj_.setModelCoefficients (table_coefficients_);
	proj_.filter (table_projected);
	table_projected_.reset (new pcl::PointCloud<Point> (table_projected));
	
	// ---[ Estimate the convex hull
	pcl::PointCloud<Point> table_hull;
	hull_.setDimension(2);
	hull_.setInputCloud (table_projected_);
	hull_.reconstruct (table_hull);
	table_hull_.reset (new pcl::PointCloud<Point> (table_hull));

	// ---[ Get the objects on top of the table
	pcl::PointIndices cloud_object_indices;
	prism_.setInputCloud (cloud_filtered_);
	prism_.setInputPlanarHull (table_hull_);
	prism_.segment (cloud_object_indices);
	
	pcl::PointCloud<Point> cloud_objects;
	pcl::ExtractIndices<Point> extract_object_indices;
	extract_object_indices.setInputCloud (cloud_filtered_);
	extract_object_indices.setIndices (boost::make_shared<const pcl::PointIndices> (cloud_object_indices));
	extract_object_indices.filter (cloud_objects);
	cloud_objects_.reset (new pcl::PointCloud<Point> (cloud_objects));
	
	if (cloud_objects.points.size () == 0)
		return false;

	// ---[ Downsample the points
#if 0
	pcl::PointCloud<Point> cloud_objects_downsampled;
	grid_objects_.setInputCloud (cloud_objects_);
	grid_objects_.filter (cloud_objects_downsampled);
	cloud_objects_downsampled_.reset (new pcl::PointCloud<Point> (cloud_objects_downsampled));

#endif
	
	// ---[ Split the objects into Euclidean clusters
	std::vector< pcl::PointIndices > object_clusters;
	cluster_.setInputCloud (cloud_objects_);
	cluster_.extract (object_clusters);
	
	
	for (size_t i = 0; i < object_clusters.size (); ++i) {
		std::vector <PointType, Eigen::aligned_allocator<PointType>> object_points;
		
		for (size_t k = 0; k < object_clusters[i].indices.size(); ++k) {
			int index = object_clusters[i].indices[k];
			Point p = cloud_objects_->points[index];
			object_points.push_back(p);
		}

		float min_dist_to_plane = FLT_MAX;
		for (int j = 0; j < object_points.size(); ++j) {
			pcl::PointXYZ pobj;
			pobj.x = object_points[j].x;
			pobj.y = object_points[j].y;
			pobj.z = object_points[j].z;

			min_dist_to_plane = std::min(plane().distanceToPlane(pobj), min_dist_to_plane);
		}

		if (min_dist_to_plane > m_max_dist_to_plane)
			continue;

		m_object_clusters.push_back(object_points);
	}
	
	return true;
}

} /* ope */