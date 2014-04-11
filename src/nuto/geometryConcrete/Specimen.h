/*
 * Specimen.h
 *
 *  Created on: 9 Apr 2014
 *      Author: ttitsche
 */

#ifndef SPECIMEN_H_
#define SPECIMEN_H_

#include "nuto/math/FullVector.h"

namespace NuTo
{
//! @brief ... class for Specimen
class Specimen
{
public:

	//! @brief ... constructor
	//! @param rBoundingBox ... bounding box of the specimen
	//! @param rTypeOfSpecimen ... element of enum Specimen::Type
	//! @param rIs2D ... true for 2D specimen
	Specimen(
			NuTo::FullMatrix<double,Eigen::Dynamic, Eigen::Dynamic> rBoundingBox,
			const int rTypeOfSpecimen,
			const bool rIs2D = false);

	//! @brief ... copy constructor
	Specimen(const NuTo::Specimen& rOther);

	//! @brief ... instead of hard coded 0,1,2...
	enum Type {Box=0, Dogbone=1, Cylinder=2};


	const bool Is2D() const;
	const bool IsBox() const;

	//! @brief ... getter for mBoundingBox
	const NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetBoundingBox() const;

	//! @brief ... getter for mLength
	const NuTo::FullVector<double, Eigen::Dynamic>& GetLength() const;

	//! @brief ... getter for specific element mLength[rIndex]
	const double GetLength(const int rIndex) const;

	//! @brief ... getter of mTypeOfSpecimen
	const int GetTypeOfSpecimen() const;

	//! @brief ... calculates and returns specimen volume
	const double GetVolume() const;

private:
	NuTo::FullMatrix<double,Eigen::Dynamic, Eigen::Dynamic> mBoundingBox;
	NuTo::FullVector<double,Eigen::Dynamic> mLength;
	const int mTypeOfSpecimen;
	const bool mIs2D;

	//! @brief ... calculates the length of the bounding box
	void CalculateLength();

	//! @brief ... throws if the bounding box is falsely defined
	void CheckBoundingBox();
};

} /* namespace NuTo */
#endif /* SPECIMEN_H_ */
