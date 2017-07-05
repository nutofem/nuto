#include "metamodel/MetamodelException.h"
#include "metamodel/ZeroMeanUnitVarianceTransformation.h"

// constructor
NuTo::ZeroMeanUnitVarianceTransformation::ZeroMeanUnitVarianceTransformation(unsigned int rCoordinate):
Transformation(),
mCoordinate(rCoordinate),
mMean(0.0),
mStandardDeviation(0.0)
{
}

// copy constructor
NuTo::ZeroMeanUnitVarianceTransformation::ZeroMeanUnitVarianceTransformation(const ZeroMeanUnitVarianceTransformation &other)
{
    mCoordinate = other.mCoordinate;
    mMean = other.mMean;
    mStandardDeviation = other.mStandardDeviation;
}


// build transformation
void NuTo::ZeroMeanUnitVarianceTransformation::Build(const Eigen::MatrixXd& rCoordinates)
{
    // check input
    if (rCoordinates.cols() < 2)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] number of points must be greater than one - check the number of columns of your matrix.");
	}
    if (rCoordinates.rows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }
    
    // calculate mean
    this->mMean = 0.0;
    const double *dataPtr = &rCoordinates.data()[mCoordinate];
    for (int count=0; count<rCoordinates.cols(); count++)
	{
        this->mMean += *dataPtr;
        dataPtr+=rCoordinates.rows();
	}
    this->mMean /= rCoordinates.cols();
    
    // calculate variance
    double variance = 0.0;
    dataPtr = &rCoordinates.data()[mCoordinate];
    for (int count=0; count<rCoordinates.cols(); count++)
	{
        double delta = *dataPtr - this->mMean;
        variance += delta * delta;
        dataPtr+=rCoordinates.rows();
    }
    variance /= rCoordinates.cols() - 1;
    this->mStandardDeviation = sqrt(variance);
    if(this->mStandardDeviation < 1e-12)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] the standard deviation is almost zero");
    }
}

// transformation
void NuTo::ZeroMeanUnitVarianceTransformation::TransformForward(Eigen::MatrixXd& rCoordinates)const
{
    // check input
    if (rCoordinates.cols() == 0)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformForward] number of points must be greater than zero - check the number of columns of your matrix.");
	}
    if (rCoordinates.rows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformForward] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }

    // transform coordinates
    double *dataPtr =  &rCoordinates.data()[mCoordinate];
    for (int count=0; count<rCoordinates.cols(); count++)
	{
	    *dataPtr = (*dataPtr - this->mMean)/this->mStandardDeviation;
        dataPtr+=rCoordinates.rows();
	}
}

// back transformation
void NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward(Eigen::MatrixXd& rCoordinates)  const
{
    // check input
    if (rCoordinates.cols() == 0)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward] number of points must be greater than zero - check the number of columns of your matrix.");
	}
    if (rCoordinates.rows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }

    // transform coordinates
    double *dataPtr =  &rCoordinates.data()[mCoordinate];
    for (int count=0; count<rCoordinates.cols(); count++)
	{
	    *dataPtr = *dataPtr * this->mStandardDeviation + this->mMean;
        dataPtr+=rCoordinates.rows();
	}
}

