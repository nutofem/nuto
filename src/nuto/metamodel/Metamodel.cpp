// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/metamodel/MinMaxTransformation.h"
#include "nuto/metamodel/ZeroMeanUnitVarianceTransformation.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/metamodel/Metamodel.h"

// constructor
NuTo::Metamodel::Metamodel() : NuTo::NuToObject()
{
	// init random number generator with milliseconds from ..
	dsfmt_init_gen_rand(&mRandomNumberGenerator, time (NULL));
}

void NuTo::Metamodel::AppendMinMaxTransformationInput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimInput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationInput - coordinate is out of range (larger than dimInput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

void NuTo::Metamodel::AppendMinMaxTransformationInput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimInput(); count++)
		AppendMinMaxTransformationInput(count, rMin, rMax);
}

void NuTo::Metamodel::AppendMinMaxTransformationOutput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimOutput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

void NuTo::Metamodel::AppendMinMaxTransformationOutput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimOutput(); count++)
		AppendMinMaxTransformationOutput(count, rMin, rMax);
}

// add zero mean, unit variance transformation to inputs
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationInput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimInput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationInput(count);
    }
}
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimInput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimInput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

// add zero mean, unit variance transformation to outputs
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationOutput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimOutput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationOutput(count);
    }
}
void NuTo::Metamodel::AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimOutput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

NuTo::FullMatrix<double> NuTo::Metamodel::GetOriginalSupportPointsInput()const
{
    return mSupportPoints.GetOrigSupportPointsInput();
}

NuTo::FullMatrix<double> NuTo::Metamodel::GetOriginalSupportPointsOutput()const
{
    return mSupportPoints.GetOrigSupportPointsOutput();
}

NuTo::FullMatrix<double> NuTo::Metamodel::GetTransformedSupportPointsInput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsInput();
}

NuTo::FullMatrix<double> NuTo::Metamodel::GetTransformedSupportPointsOutput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsOutput();
}

void NuTo::Metamodel::SetSupportPoints(int rDimInput, int rDimOutput, FullMatrix<double> rInputCoordinates, FullMatrix<double> rOutputCoordinates)
{
    if (rDimInput!=rInputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of input  must be equal to number of rows in the input matrix.");

    if (rDimOutput!=rOutputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of output  must be equal to number of rows in the output matrix.");

    if (rOutputCoordinates.GetNumColumns()!=rInputCoordinates.GetNumColumns())
	    throw MetamodelException("Metamodel::SetSupportPoints - number of samples (number of columns) in input and output matrix must be identical .");
		
	mSupportPoints.SetSupportPoints(rInputCoordinates,rOutputCoordinates);
}

void NuTo::Metamodel::BuildTransformation()
{
	if (!mSupportPoints.IsTransformationBuild())
		mSupportPoints.BuildTransformation();
}

void NuTo::Metamodel::InitRandomNumberGenerator(int rSeed)
{
    dsfmt_init_gen_rand(&mRandomNumberGenerator, rSeed);  
}

double NuTo::Metamodel::RandomDouble()
{
    return dsfmt_genrand_close_open(&mRandomNumberGenerator);
}

void NuTo::Metamodel::Solve(const FullMatrix<double>& rInputCoordinates, FullMatrix<double>& rOutputCoordinates)const
{
    //apply transformation of inputs
	NuTo::FullMatrix<double> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);

    //solve the submodule
    SolveTransformed(rInputCoordinatesTransformed, rOutputCoordinates);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
}

void NuTo::Metamodel::SolveConfidenceInterval(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates,
                             NuTo::FullMatrix<double>& rOutputCoordinatesMin, NuTo::FullMatrix<double>& rOutputCoordinatesMax)const
{
    //apply transformation of inputs
	NuTo::FullMatrix<double> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);
    
    //solve the submodule
    SolveConfidenceIntervalTransformed(rInputCoordinatesTransformed, rOutputCoordinates, rOutputCoordinatesMin, rOutputCoordinatesMax);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMin);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMax);
    
}
void NuTo::Metamodel::Build()
{
    if (mSupportPoints.GetDimOutput()<1)
        throw MetamodelException("NuTo::Metamodel::Build - number of outputs must be positive - set training data first.");
    
    if (mSupportPoints.GetDimInput()<1)
        throw MetamodelException("NuTo::Metamodel::Build - number of inputs must be positive - set training data first.");

    if (!mSupportPoints.IsTransformationBuild())
        mSupportPoints.BuildTransformation();
        
    BuildDerived();
}

void NuTo::Metamodel::Info()const
{
	mSupportPoints.Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Metamodel::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Metamodel::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Metamodel::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize Metamodel" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(NuToObject)
       & BOOST_SERIALIZATION_NVP(mSupportPoints);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize Metamodel" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Metamodel)
#endif  // ENABLE_SERIALIZATION
