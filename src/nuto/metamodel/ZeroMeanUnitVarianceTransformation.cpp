// $Id$

#include "nuto/metamodel/MetamodelException.h"
#include "nuto/metamodel/ZeroMeanUnitVarianceTransformation.h"

// constructor
NuTo::ZeroMeanUnitVarianceTransformation::ZeroMeanUnitVarianceTransformation(unsigned int rCoordinate) : Transformation()
{
    if(rCoordinate < 0)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::ZeroMeanUnitVarianceTransformation] coordinate must be a positive value.");
    }
    this->mCoordinate = rCoordinate; 
    this->mMean = 0.0;
    this->mStandardDeviation = 0.0;
}

// build transformation
void NuTo::ZeroMeanUnitVarianceTransformation::Build(const FullMatrix<double>& rCoordinates)
{
    // check input
    if (rCoordinates.GetNumColumns() < 2)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] number of points must be greater than one - check the number of columns of your matrix.");
	}
    if (rCoordinates.GetNumRows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }
    
    // calculate mean
    this->mMean = 0.0;
    const double *dataPtr = &rCoordinates.mEigenMatrix.data()[mCoordinate];
    for (int count=0; count<rCoordinates.GetNumColumns(); count++)
	{
        this->mMean += *dataPtr;
        dataPtr+=rCoordinates.GetNumRows();
	}
    this->mMean /= rCoordinates.GetNumColumns();
    
    // calculate variance
    double variance = 0.0;
    dataPtr = &rCoordinates.mEigenMatrix.data()[mCoordinate];
    for (int count=0; count<rCoordinates.GetNumColumns(); count++)
	{
        double delta = *dataPtr - this->mMean;
        variance += delta * delta;
        dataPtr+=rCoordinates.GetNumRows();
    }
    variance /= rCoordinates.GetNumColumns() - 1;
    this->mStandardDeviation = sqrt(variance);
    if(this->mStandardDeviation < 1e-12)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::Build] the standard deviation is almost zero");
    }
}

// transformation
void NuTo::ZeroMeanUnitVarianceTransformation::TransformForward(FullMatrix<double>& rCoordinates)const
{
    // check input
    if (rCoordinates.GetNumColumns() == 0)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformForward] number of points must be greater than zero - check the number of columns of your matrix.");
	}
    if (rCoordinates.GetNumRows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformForward] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }

    // transform coordinates
    double *dataPtr =  &rCoordinates.mEigenMatrix.data()[mCoordinate];
    for (int count=0; count<rCoordinates.GetNumColumns(); count++)
	{
	    *dataPtr = (*dataPtr - this->mMean)/this->mStandardDeviation;
        dataPtr+=rCoordinates.GetNumRows();
	}
}

// back transformation
void NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward(FullMatrix<double>& rCoordinates)  const
{
    // check input
    if (rCoordinates.GetNumColumns() == 0)
	{
	    throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward] number of points must be greater than zero - check the number of columns of your matrix.");
	}
    if (rCoordinates.GetNumRows() <= this->mCoordinate)
    {
        throw MetamodelException("[NuTo::ZeroMeanUnitVarianceTransformation::TransformBackward] coordinate to be transformed is out of range - check the number of rows of your Matrix.");
    }

    // transform coordinates
    double *dataPtr =  &rCoordinates.mEigenMatrix.data()[mCoordinate];
    for (int count=0; count<rCoordinates.GetNumColumns(); count++)
	{
	    *dataPtr = *dataPtr * this->mStandardDeviation + this->mMean;
        dataPtr+=rCoordinates.GetNumRows();
	}
}