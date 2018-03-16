/*
 * Specimen.h
 *
 *  Created on: 9 Apr 2014
 *      Author: ttitsche
 */

#pragma once

#include <Eigen/Core>

namespace NuTo
{

//! @brief class for Specimen
class Specimen
{
public:
    enum eSpecimenType
    {
        Box = 0,
        Dogbone = 1,
        Cylinder = 2
    };

    //! @brief constructor
    //! @param rBoundingBox bounding box of the specimen
    //! @param rTypeOfSpecimen element of enum Specimen::Type
    Specimen(Eigen::MatrixXd rBoundingBox, const eSpecimenType rTypeOfSpecimen);

    //! @brief copy constructor
    Specimen(const NuTo::Specimen& rOther);

    bool IsBox() const;

    //! @brief getter for mBoundingBox
    const Eigen::MatrixXd& GetBoundingBox() const;

    //! @brief getter for mLength
    const Eigen::VectorXd& GetLength() const;

    //! @brief getter for specific element mLength[rIndex]
    double GetLength(const int rIndex) const;

    //! @brief getter of mTypeOfSpecimen
    eSpecimenType GetTypeOfSpecimen() const;

    //! @brief calculates and returns specimen volume
    double GetVolume() const;

private:
    Eigen::MatrixXd mBoundingBox;
    Eigen::VectorXd mLength;
    const eSpecimenType mTypeOfSpecimen;

    //! @brief calculates the length of the bounding box
    void CalculateLength();

    //! @brief throws if the bounding box is falsely defined
    void CheckBoundingBox();
};

} /* namespace NuTo */
