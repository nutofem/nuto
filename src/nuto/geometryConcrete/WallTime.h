/*
 * WallTime.h
 *
 *  Created on: 17 Mar 2014
 *      Author: ttitsche
 */

#pragma once

#include <string>


namespace NuTo
{

class WallTime
{
public:
	/**
	 *
	 * @return
	 */
	static double Get();

	//! @brief ... returns a TimeStamp in YYYY-MM-DD_HH:MM:SS format
	static std::string TimeStamp();
};

} /* namespace NuTo */
