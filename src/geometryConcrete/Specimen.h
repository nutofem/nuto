/*
 * Specimen.h
 *
 *  Created on: 9 Apr 2014
 *      Author: ttitsche
 */

#pragma once

#include <eigen3/Eigen/Core>
#include "math/FullMatrix.h"

namespace NuTo
{

//! @brief ... class for Specimen
class Specimen
{
public:

    enum eSpecimenType
    {
        Box=0,
        Dogbone=1,
        Cylinder=2
    };

	//! @brief ... constructor
	//! @param rBoundingBox ... bounding box of the specimen
	//! @param rTypeOfSpecimen ... element of enum Specimen::Type
	Specimen(
			NuTo::FullMatrix<double,Eigen::Dynamic, Eigen::Dynamic> rBoundingBox,
            const eSpecimenType rTypeOfSpecimen);

	//! @brief ... copy constructor
	Specimen(const NuTo::Specimen& rOther);

	const bool IsBox() const;

	//! @brief ... getter for mBoundingBox
	const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetBoundingBox() const;

	//! @brief ... getter for mLength
	const Eigen::VectorXd& GetLength() const;

	//! @brief ... getter for specific element mLength[rIndex]
	const double GetLength(const int rIndex) const;

	//! @brief ... getter of mTypeOfSpecimen
	const eSpecimenType GetTypeOfSpecimen() const;

	//! @brief ... calculates and returns specimen volume
	const double GetVolume() const;

private:
	NuTo::FullMatrix<double,Eigen::Dynamic, Eigen::Dynamic> mBoundingBox;
	Eigen::VectorXd mLength;
	const eSpecimenType mTypeOfSpecimen;

	//! @brief ... calculates the length of the bounding box
	void CalculateLength();

	//! @brief ... throws if the bounding box is falsely defined
	void CheckBoundingBox();
};

} /* namespace NuTo */
