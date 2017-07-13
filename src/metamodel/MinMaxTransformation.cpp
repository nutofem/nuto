#include "base/Exception.h"
#include "metamodel/MinMaxTransformation.h"


// constructor
NuTo::MinMaxTransformation::MinMaxTransformation(unsigned int rCoordinate, double rLb, double rUb)
    : NuTo::Transformation()
    , mCoordinate(rCoordinate)
    , mMin(0.0)
    , mMax(0.0)
    , mUb(rUb)
    , mLb(rLb)
{
}

NuTo::MinMaxTransformation::MinMaxTransformation(const MinMaxTransformation& rOther)
{
    mCoordinate = rOther.mCoordinate;
    mLb = rOther.mLb;
    mUb = rOther.mUb;
    mMin = rOther.mMin;
    mMax = rOther.mMax;
}


void NuTo::MinMaxTransformation::Build(const Eigen::MatrixXd& rCoordinates)
{
    if (rCoordinates.cols() == 0)
    {
        throw Exception("MinMaxTransformation::build - numberOfPoints must be greater than zero");
    }
    if (rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::build - coordinate to be transformed is out of range - check "
                                 "the dimension of your Matrix.");
    }
    const double* theptr = &rCoordinates.data()[mCoordinate];
    mMin = *theptr;
    mMax = *theptr;
    for (int count = 1; count < rCoordinates.cols(); count++)
    {
        theptr += rCoordinates.rows();
        if (*theptr < mMin)
            mMin = *theptr;
        else
        {
            if (*theptr > mMax)
                mMax = *theptr;
        }
    }
}

void NuTo::MinMaxTransformation::TransformForward(Eigen::MatrixXd& rCoordinates) const
{
    if (rCoordinates.cols() == 0)
    {
        throw Exception("MinMaxTransformation::TransformForward - numberOfPoints must be greater than zero");
    }
    if (rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::TransformForward - coordinate to be transformed is out of "
                                 "range - check the dimension of your Matrix.");
    }
    double* theptr = &rCoordinates.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax - mMin;

    if (deltaBound == 0)
    {
        throw Exception("MinMaxTransformation::TransformForward - delta of prescribed bounds equal to zero");
    }

    if (deltaValue == 0)
    {
        throw Exception("MinMaxTransformation::TransformForward - interval between min and max value of given "
                                 "points has size zero");
    }

    for (int count = 0; count < rCoordinates.cols(); count++, theptr += rCoordinates.rows())
    {
        *theptr = mLb + (*theptr - mMin) / deltaValue * deltaBound;
    }
}

void NuTo::MinMaxTransformation::TransformBackward(Eigen::MatrixXd& rCoordinates) const
{
    if (rCoordinates.cols() == 0)
    {
        throw Exception("MinMaxTransformation::TransformBackward - numberOfPoints must be greater than zero");
    }
    if (rCoordinates.rows() <= mCoordinate)
    {
        throw Exception("MinMaxTransformation::TransformBackward - coordinate to be transformed is out of "
                                 "range - check the dimension of your Matrix.");
    }
    double* theptr = &rCoordinates.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax - mMin;

    if (deltaBound == 0)
    {
        throw Exception("MinMaxTransformation::TransformBackward - delta of prescribed bounds equal to zero");
    }

    if (deltaValue == 0)
    {
        throw Exception("MinMaxTransformation::TransformBackward - interval between min and max value of "
                                 "given points has size zero");
    }

    for (int count = 0; count < rCoordinates.cols(); count++, theptr += rCoordinates.rows())
    {
        *theptr = mMin + (*theptr - mLb) / deltaBound * deltaValue;
    }
}
