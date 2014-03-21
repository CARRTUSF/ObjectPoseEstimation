/*
 *  Object Pose Estimation (OPE) - www.cse.usf.edu/kkduncan/ope
 *  Copyright (c) 2013, Kester Duncan
 *
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <boost/math/special_functions/factorials.hpp>

#include "SQTypes.h"
#include "InertiaCalculations.h"


namespace ope {


double InertiaCalculations::betaFunction (const double& x, const double& y) {
	double result = 0.0;

	result = (boost::math::factorial<double>(static_cast<unsigned>(x - 1)) * 
		boost::math::factorial<double>(static_cast<unsigned>(y - 1))) / 
		boost::math::factorial<double>(static_cast<unsigned>(x + y - 1));

	return result;
}


double InertiaCalculations::xMomentOfInertia(SQParameters& sq) {
	double result = 0.0;

	result = 0.5 * sq.a1 * sq.a2 * sq.a3 * sq.e1 * sq.e2 * ((sq.a2 * sq.a2) * betaFunction(1.5 * sq.e2, 0.5 * sq.e2) *
		betaFunction(0.5 * sq.e1, (2 * sq.e1 + 1)) + (4 * sq.a3 * sq.a3) * betaFunction(0.5 * sq.e2, (0.5 * sq.e2 + 1)) *
		betaFunction(1.5 * sq.e1, sq.e1 + 1));

	return result;
}


double InertiaCalculations::yMomentOfInertia(SQParameters& sq) {
	double result = 0.0;

	result = 0.5 * sq.a1 * sq.a2 * sq.a3 * sq.e1 * sq.e2 * ((sq.a1 * sq.a1) * betaFunction(1.5 * sq.e2, 0.5 * sq.e2) *
		betaFunction(0.5 * sq.e1, (2 * sq.e1 + 1)) + (4 * sq.a3 * sq.a3) * betaFunction(0.5 * sq.e2, (0.5 * sq.e2 + 1)) *
		betaFunction(1.5 * sq.e1, sq.e1 + 1));

	return result;
}


double InertiaCalculations::zMomentOfInertia(SQParameters& sq) {
	double result = 0.0;

	result = 0.5 * sq.a1 * sq.a2 * sq.a3 * sq.e1 * sq.e2 * ((sq.a1 * sq.a1) + (sq.a2 * sq.a2)) * 
		betaFunction(1.5 * sq.e2, 0.5 * sq.e2) * betaFunction(0.5 * sq.e1, (2 * sq.e1 + 1));

	return result;
}


void InertiaCalculations::calculateNewPrincipalAxis (SQParameters& sq) {
	double xAxisInertia = xMomentOfInertia(sq);
	double yAxisInertia = yMomentOfInertia(sq);
	double zAxisInertia = zMomentOfInertia(sq);
	double minInertiaValue = 99999.0;
	int minInertiaAxis = -1;

	if (xAxisInertia < minInertiaValue) {
		minInertiaValue = xAxisInertia;
		minInertiaAxis = 0;
	}

	if (yAxisInertia < minInertiaValue) {
		minInertiaValue = yAxisInertia;
		minInertiaAxis = 1;
	}

	if (zAxisInertia < minInertiaValue) {
		minInertiaValue = zAxisInertia;
		minInertiaAxis = 2;
	}

	if (minInertiaAxis == -1) {
		minInertiaAxis = sq.principalAxis;
	}

	sq.principalAxis = minInertiaAxis;
}


void InertiaCalculations::calculateMajorAxis (SQParameters& sq) {
	double xAxisInertia = xMomentOfInertia(sq);
	double yAxisInertia = yMomentOfInertia(sq);
	double zAxisInertia = zMomentOfInertia(sq);
	double minInertiaValue = 99999.0;
	int minInertiaAxis = -1;

	if (xAxisInertia < minInertiaValue) {
		minInertiaValue = xAxisInertia;
		minInertiaAxis = 0;
	}

	if (yAxisInertia < minInertiaValue) {
		minInertiaValue = yAxisInertia;
		minInertiaAxis = 1;
	}

	if (zAxisInertia < minInertiaValue) {
		minInertiaValue = zAxisInertia;
		minInertiaAxis = 2;
	}

	if (minInertiaAxis == -1) {
		minInertiaAxis = sq.principalAxis;
	}

	sq.majorAxis = minInertiaAxis;
}


void InertiaCalculations::calculateMinorAxis (SQParameters& sq) {
	double xAxisInertia = xMomentOfInertia(sq);
	double yAxisInertia = yMomentOfInertia(sq);
	double zAxisInertia = zMomentOfInertia(sq);
	double maxInertiaValue = -99999.0;
	int maxInertiaAxis = -1;

	if (xAxisInertia > maxInertiaValue) {
		maxInertiaValue = xAxisInertia;
		maxInertiaAxis = 0;
	}

	if (yAxisInertia > maxInertiaValue) {
		maxInertiaValue = yAxisInertia;
		maxInertiaAxis = 1;
	}

	if (zAxisInertia > maxInertiaValue) {
		maxInertiaValue = zAxisInertia;
		maxInertiaAxis = 2;
	}

	if (maxInertiaAxis == -1) {
		maxInertiaAxis = sq.principalAxis;
	}

	sq.minorAxis = maxInertiaAxis;

}


} /* ope namespace */

