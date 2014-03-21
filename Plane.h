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
 *	\file Plane.h
 *	\author Kester Duncan
 *  \note Adapted from ntk's <code>Plane</code> object
 */
#pragma once
#ifndef PLANE_H__
#define PLANE_H__


// Forward declaration
namespace pcl {
	class PointXYZ;

}


namespace ope
{


/**
 * \brief Defines the properties of a plane
 */
class Plane {
	/// Plane parameters
	double a, b, c, d;

public:
	Plane(double a, double b, double c, double d)
		: a(a), b(b), c(c), d(d)
	{}

	/**
	 * \brief Construct a plane from a normal vector and a point.
	 */
	Plane(const pcl::PointXYZ& normal, const pcl::PointXYZ& p);

	/**
	 * \brief Default Constructor
	 */
	Plane() : a(0), b(0), c(0), d(0)
	{}

	/**
	 * \brief Determines whether or not this is a valid plane
	 * \return true if it is valid
	 */
	bool isValid() const;

	/**
	 * \brief Gets the normal vector that defines this plane
	 * \returns the x, y, & z values of the plane's normal
	 */
	pcl::PointXYZ normal() const;

	/**
	 * \brief Sets the parameters of the plane
	 */
	void set (double a_, double b_, double c_, double d_) { 
		a = a_; 
		b = b_; 
		c = c_; 
		d = d_; 
	}

	/**
	 * \brief Determines whether this plane intersects with the specified line given by the parameters
	 * \param <p1> the first point that defines the line
	 * \param <p2> the second point that defines the line
	 */
	pcl::PointXYZ intersectionWithLine (const pcl::PointXYZ& p1, const pcl::PointXYZ& p2) const;


	/**
	 * \brief Determines a point's distance from the plane
	 * \return Distance from the plane
	 */
	float distanceToPlane(const pcl::PointXYZ& p) const;


};


} // ope

#endif /* PLANE_H__ */
