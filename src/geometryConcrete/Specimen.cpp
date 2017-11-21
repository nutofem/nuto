/*
 * Specimen.cpp
 *
 *  Created on: 9 Apr 2014
 *      Author: ttitsche
 */

#include "base/Exception.h"
#include "geometryConcrete/Specimen.h"

NuTo::Specimen::Specimen(Eigen::MatrixXd rBoundingBox, const eSpecimenType rTypeOfSpecimen)
    : mBoundingBox(rBoundingBox)
    , mTypeOfSpecimen(rTypeOfSpecimen)
{
    CheckBoundingBox();
    CalculateLength();
}

NuTo::Specimen::Specimen(const Specimen& rOther)
    : mBoundingBox(rOther.mBoundingBox)
    , mTypeOfSpecimen(rOther.mTypeOfSpecimen)
{
    CalculateLength();
}

double NuTo::Specimen::GetVolume() const
{
    if (mTypeOfSpecimen != 0 && mTypeOfSpecimen != 2)
        throw Exception(__PRETTY_FUNCTION__, "specimen type not implemented.");

    double Vspecimen;

    switch (mTypeOfSpecimen)
    {
    case 0:
        // box
        Vspecimen = mLength[0] * mLength[1] * mLength[2];
        break;
    case 1:
    {
        Vspecimen = mLength[0] * mLength[1] * mLength[2];
        double D = mLength[0];
        if (std::abs(static_cast<double>(mLength[1] - 1.5 * D)) > 1e-10)
            throw Exception(__PRETTY_FUNCTION__,
                            "for the dog bone specimen, the y dimension should be 1.5 times the x dimension.");
        // subtract the circles
        double radius = 0.725 * D;
        double deltaAngle = 2. * 0.2 / 0.525;
        Vspecimen -= 2. * (deltaAngle / (2. * M_PI) * M_PI * radius * radius) * mLength[2];
    }
    break;
    case 2:
    {
        if (std::abs(static_cast<double>(mLength[0] - mLength[1])) > 1e-10)
            throw Exception(__PRETTY_FUNCTION__,
                            "for the cylindern, the x and y dimension should be identical (Diameter).");
        double D = mLength[0];
        if (D < 1e-10)
            throw Exception(__PRETTY_FUNCTION__, "for the cylindern, the x,y dimension should be positive (Diameter).");

        Vspecimen = M_PI * 0.25 * D * D * mLength[2];
    }
    break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "specimen type not implemented.");
    }
    if (Vspecimen < 1.0e-14)
    {
        throw Exception(__PRETTY_FUNCTION__, "negative volume of the box.");
    }
    return Vspecimen;
}

void NuTo::Specimen::CalculateLength()
{
    mLength.resize(3);

    int dimensionsToCheck = 3;

    for (int i = 0; i < dimensionsToCheck; i++)
    {
        mLength[i] = mBoundingBox(i, 1) - mBoundingBox(i, 0);
        if (mLength[i] <= 0)
            throw Exception(__PRETTY_FUNCTION__, "box dimensions should be not negative.");
    }
}

const Eigen::MatrixXd& NuTo::Specimen::GetBoundingBox() const
{
    return mBoundingBox;
}

const Eigen::VectorXd& NuTo::Specimen::GetLength() const
{
    return mLength;
}

NuTo::Specimen::eSpecimenType NuTo::Specimen::GetTypeOfSpecimen() const
{
    return mTypeOfSpecimen;
}

bool NuTo::Specimen::IsBox() const
{
    return mTypeOfSpecimen == eSpecimenType::Box;
}

double NuTo::Specimen::GetLength(const int rIndex) const
{
    return mLength[rIndex];
}

void NuTo::Specimen::CheckBoundingBox()
{
    if (mBoundingBox.rows() != 3 && mBoundingBox.cols() != 2)
        throw Exception(__PRETTY_FUNCTION__, "bounding box has to have the dimension [3,2]");

    //  if (rBoundingBox.GetValue(0, 0) != 0. || rBoundingBox.GetValue(1, 0) != 0. || rBoundingBox.GetValue(2, 0) != 0.)
    //      throw Exception("[NuTo::ParticleCreator::CheckBoundingBox] bounding box has to start at (0.,0.,0.)");
}
