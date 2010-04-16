/*******************************************************************************
Bauhaus-Universitaet Weimar
Author: Joerg F. Unger,  Septermber 2009
*******************************************************************************/

#include "nuto/metamodel/MetamodelException.h"
#include "nuto/metamodel/MinMaxTransformation.h"

using namespace NuTo;

void MinMaxTransformation::Build(const FullMatrix<double>& rCoordinates)
{
    if ( rCoordinates.GetNumColumns() == 0)
	{
	    throw MetamodelException("MinMaxTransformation::build - numberOfPoints must be greater than zero");
	}
    if ( rCoordinates.GetNumRows() <= mCoordinate)
    {
        throw MetamodelException("MinMaxTransformation::build - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    const double *theptr = &rCoordinates.mEigenMatrix.data()[mCoordinate];
	mMin = *theptr;
	mMax = *theptr;
    for (int count=1; count<rCoordinates.GetNumColumns(); count++)
	{
        theptr+=rCoordinates.GetNumRows();
		if (*theptr<mMin)
			mMin = *theptr;
		else
		{
		    if (*theptr>mMax)
			    mMax = *theptr;
		}
	}
}

void MinMaxTransformation::TransformForward(FullMatrix<double>& rCoordinates)const
{
    if ( rCoordinates.GetNumColumns() == 0)
    {
        throw MetamodelException("MinMaxTransformation::TransformForward - numberOfPoints must be greater than zero");
    }
    if ( rCoordinates.GetNumRows() <= mCoordinate)
    {
        throw MetamodelException("MinMaxTransformation::TransformForward - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    double *theptr =  &rCoordinates.mEigenMatrix.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax-mMin;
	
	if (deltaBound==0 )
	{
        throw MetamodelException("MinMaxTransformation::TransformForward - delta of prescribed bounds equal to zero");
	}

	if (deltaValue==0)
	{
        throw MetamodelException("MinMaxTransformation::TransformForward - interval between min and max value of given points has size zero");
	}

    for (int count=0; count<rCoordinates.GetNumColumns(); count++,theptr+=rCoordinates.GetNumRows())
	{
	    *theptr = mLb + (*theptr-mMin)/deltaValue*deltaBound;
	}
}

void MinMaxTransformation::TransformBackward(FullMatrix<double>& rCoordinates)  const
{
    if ( rCoordinates.GetNumColumns() == 0)
    {
        throw MetamodelException("MinMaxTransformation::TransformBackward - numberOfPoints must be greater than zero");
    }
    if ( rCoordinates.GetNumRows() <= mCoordinate)
    {
        throw MetamodelException("MinMaxTransformation::TransformBackward - coordinate to be transformed is out of range - check the dimension of your Matrix.");
    }
    double *theptr =  &rCoordinates.mEigenMatrix.data()[mCoordinate];
    double deltaBound = mUb - mLb;
    double deltaValue = mMax-mMin;
	
	if (deltaBound==0 )
	{
        throw MetamodelException("MinMaxTransformation::TransformBackward - delta of prescribed bounds equal to zero");
	}

	if (deltaValue==0)
	{
        throw MetamodelException("MinMaxTransformation::TransformBackward - interval between min and max value of given points has size zero");
	}

    for (int count=0; count<rCoordinates.GetNumColumns(); count++,theptr+=rCoordinates.GetNumRows())
	{
	    *theptr = mMin + (*theptr-mLb)/deltaBound*deltaValue;
	}
}
