/*
 * Point.cpp
 *
 *  Created on: Jun 22, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <cmath>
#include <limits>

#include "Point.h"

bool operator<=(const GEOLIB::Point& p0, const GEOLIB::Point& p1)
{
	double tol(sqrt(std::numeric_limits<double>::min()));

	if (fabs(p0[0] - p1[0]) > tol * fabs(p0[0]))
	{
		if (p0[0] < p1[0])
			return true;
		else
			return false;
	}
	else
	{
		// assume p0[0] == p1[0]
		if (fabs(p0[1] - p1[1]) > tol * fabs(p0[1]))
		{
			if (p0[1] < p1[1])
				return true;
			else
				return false;
		}
		else
		{
			// assume p0[1] == p1[1] and p0[0] == p1[0]
			if (p0[2] < p1[2])
				return true;
			else
				return false;
		}
	}
}

namespace GEOLIB
{
bool lessX(GEOLIB::Point const& p0, GEOLIB::Point const& p1)
{
	if (p0[0] <= p1[0]) return true;
	return false;
}

bool lessY(GEOLIB::Point const& p0, GEOLIB::Point const& p1)
{
	if (p0[1] <= p1[1]) return true;
	return false;
}

bool lessZ(GEOLIB::Point const& p0, GEOLIB::Point const& p1)
{
	if (p0[2] <= p1[2]) return true;
	return false;
}
}  // end namespace GEOLIB
