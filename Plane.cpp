#include <pcl/point_types.h>
#include "Plane.h"

using namespace ope;


/// Normalize an XYZ vector
static void normalize(pcl::PointXYZ& v) {
	float denom = 1.0f / std::sqrtf((v.x * v.x) + (v.y * v.y) + (v.z * v.z));

	v.x *= denom;
	v.y *= denom;
	v.z *= denom;
}


Plane :: Plane(const pcl::PointXYZ& normal, const pcl::PointXYZ& p) {
	a = normal.x;
	b = normal.y;
	c = normal.z;
	d = -(normal.x * p.x + normal.y * p.y + normal.z *p.z);
}

pcl::PointXYZ Plane :: normal() const {
	pcl::PointXYZ n(a,b,c);
	normalize(n);
	return n;
}

bool Plane :: isValid() const {
	const float eps = 1e-15;
	return (std::abs(a) > eps) || (std::abs(b) > eps) || (std::abs(c) > eps);
}


pcl::PointXYZ Plane :: intersectionWithLine (const pcl::PointXYZ& p1, const pcl::PointXYZ& p2) const {
	double u = a * p1.x + b * p1.y + c * p1.z + d;
	u /= a * (p1.x - p2.x) + b * (p1.y - p2.y) + c * (p1.z - p2.z);

	pcl::PointXYZ r;  
	r.x = p1.x + u * (p2.x - p1.x);
	r.y = p1.y + u * (p2.y - p1.y);
	r.z = p1.z + u * (p2.z - p1.z);

	return r;
}


float Plane :: distanceToPlane(const pcl::PointXYZ& p) const {
	float v = a * p.x + b * p.y + c * p.z + d;
	v /= sqrt(a * a + b * b + c * c);
	return std::abs(v);
}