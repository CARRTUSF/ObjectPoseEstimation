ObjectPoseEstimation
====================

Estimates the position and orientation of an object represented as a 3D point cloud using superquadrics. 
This is the implementation of the algorithm detailed in the ICRA 2013 submission by K. Duncan et al. entitled 
"Multi-scale Superquadric Fitting for Efficient Shape and Pose Recovery of Unknown Objects." The article in question 
can be found here --> http://www.cse.usf.edu/~kkduncan/research/DuncanICRA2013.pdf


USAGE
=====

See OPEMain.cpp for a demonstration of how to use this code. There are two main functions to use:

  a) calculateObjectPose(pcl::PointCloud\<pcl::PointXYZRGB\>& selectedObjectPtCloud)
  
    - Found in ObjectPoseEstimator.[h,cpp]. This function is responsible for calculating the pose
      of a segmented object represented by the PCL point cloud provided. This function assumes that
      the cloud follows the Kinect optical frame. This is then transformed into world coordinates
      before superquadric processing. The function returns an instance of the SQParameters class
      which possesses the position and orientation of the object. See the comments in SQParameters.h
      for an explanation.
      
  b) run()
  
    - Also found ObjectPoseEstimator.[h,cpp]. This function assumes that there is a connected Kinect that
      faces a table-top of objects within viewing range and that the OpenNI drivers are installed. It 
      executes as follows: 
      
        1) Capture a point cloud from the Kinect.
        2) Extract the table-top objects from the captured point cloud.
        3) Present a viewing window to the user in order for them to see the index of the object they
           intend to choose.
        4) Asks the user for the index of the object they wish to choose.
        5) Performs pose estimation on the object point cloud the user indicated.
        6) Return the SQParameters instance.
        
        
NOTE
====

Regular updates to this code would be carried out in the upcoming weeks, so please keep an eye out. There
may be a few bugs in this code as error conditions are not handled as yet. Therefore use at your own
risk. This code is intended for use with Visual Studio 2010 and later. However, this code can be adapted for 
any system as long as the following libraries are available: Boost, PCL and all libraries required by PCL (Eigen, VTK, Flann, Qt, QHull, OpenNI), and OpenCV.

For any questions or comments, please contact kkduncan@cse.usf.edu
        
