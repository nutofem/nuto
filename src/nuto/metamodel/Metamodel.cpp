#include "nuto/metamodel/Metamodel.h"
#include "nuto/metamodel/MinMaxTransformation.h"
#include "nuto/metamodel/ZeroMeanUnitVarianceTransformation.h"
#include "nuto/math/FullMatrix.h"

using namespace NuTo;

void Metamodel::AppendMinMaxTransformationInput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimInput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationInput - coordinate is out of range (larger than dimInput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

void Metamodel::AppendMinMaxTransformationInput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimInput(); count++)
		AppendMinMaxTransformationInput(count, rMin, rMax);
}

void Metamodel::AppendMinMaxTransformationOutput(int rCoordinate, double rMin, double rMax)
{
    if (rCoordinate>=mSupportPoints.GetDimOutput())
	{
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
	    
	}
	MinMaxTransformation *newTransformation = new MinMaxTransformation(rCoordinate, rMin, rMax);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

void Metamodel::AppendMinMaxTransformationOutput(double rMin, double rMax)
{
    for (int count=0; count<mSupportPoints.GetDimOutput(); count++)
		AppendMinMaxTransformationOutput(count, rMin, rMax);
}

// add zero mean, unit variance transformation to inputs
void Metamodel::AppendZeroMeanUnitVarianceTransformationInput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimInput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationInput(count);
    }
}
void Metamodel::AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimInput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimInput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationInput(newTransformation);
}

// add zero mean, unit variance transformation to outputs
void Metamodel::AppendZeroMeanUnitVarianceTransformationOutput()
{
    for (int count = 0; count < this->mSupportPoints.GetDimOutput(); count++)
    {
		this->AppendZeroMeanUnitVarianceTransformationOutput(count);
    }
}
void Metamodel::AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate)
{
    if( (rCoordinate < 0) || (rCoordinate >= this->mSupportPoints.GetDimOutput()) )
    {
	    throw MetamodelException("Metamodel::AppendMinMaxTransformationOutput - coordinate is out of range (larger than dimOutput).");
    }
	ZeroMeanUnitVarianceTransformation *newTransformation = new ZeroMeanUnitVarianceTransformation(rCoordinate);
	mSupportPoints.AppendTransformationOutput(newTransformation);
}

FullMatrix<double> Metamodel::GetOriginalSupportPointsInput()const
{
    return mSupportPoints.GetOrigSupportPointsInput();
}

FullMatrix<double> Metamodel::GetOriginalSupportPointsOutput()const
{
    return mSupportPoints.GetOrigSupportPointsOutput();
}

FullMatrix<double> Metamodel::GetTransformedSupportPointsInput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsInput();
}

FullMatrix<double> Metamodel::GetTransformedSupportPointsOutput()const
{
	if (!mSupportPoints.IsTransformationBuild())
	    throw MetamodelException("Metamodel::GetTransformedSupportPoints - build the transformation first.");
    
	return mSupportPoints.GetTransformedSupportPointsOutput();
}

void Metamodel::SetSupportPoints(int rDimInput, int rDimOutput, FullMatrix<double> rInputCoordinates, FullMatrix<double> rOutputCoordinates)
{
    if (rDimInput!=rInputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of input  must be equal to number of rows in the input matrix.");

    if (rDimOutput!=rOutputCoordinates.GetNumRows())
	    throw MetamodelException("Metamodel::SetSupportPoints - dimension of output  must be equal to number of rows in the output matrix.");

    if (rOutputCoordinates.GetNumColumns()!=rInputCoordinates.GetNumColumns())
	    throw MetamodelException("Metamodel::SetSupportPoints - number of samples (number of columns) in input and output matrix must be identical .");
		
	mSupportPoints.SetSupportPoints(rInputCoordinates,rOutputCoordinates);
}

void Metamodel::BuildTransformation()
{
	if (!mSupportPoints.IsTransformationBuild())
		mSupportPoints.BuildTransformation();
}

void Metamodel::InitRandomNumberGenerator(int rSeed)
{
    dsfmt_init_gen_rand(&mRandomNumberGenerator, rSeed);  
}

double Metamodel::RandomDouble()
{
    return dsfmt_genrand_close_open(&mRandomNumberGenerator);
}

void Metamodel::Solve(const FullMatrix<double>& rInputCoordinates, FullMatrix<double>& rOutputCoordinates)const
{
    //apply transformation of inputs
    FullMatrix<double> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);

    //solve the submodule
    SolveTransformed(rInputCoordinatesTransformed, rOutputCoordinates);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
}

void Metamodel::SolveConfidenceInterval(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates, 
                             NuTo::FullMatrix<double>& rOutputCoordinatesMin, NuTo::FullMatrix<double>& rOutputCoordinatesMax)const
{
    //apply transformation of inputs
    FullMatrix<double> rInputCoordinatesTransformed = rInputCoordinates;
    mSupportPoints.TransformForwardInput(rInputCoordinatesTransformed);
    
    //solve the submodule
    SolveConfidenceIntervalTransformed(rInputCoordinatesTransformed, rOutputCoordinates, rOutputCoordinatesMin, rOutputCoordinatesMax);
    
    //apply transformation of outputs
    mSupportPoints.TransformForwardOutput(rOutputCoordinates);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMin);
    mSupportPoints.TransformForwardOutput(rOutputCoordinatesMax);
    
}
void Metamodel::Build()
{
    if (mSupportPoints.GetDimOutput()<=0)
        throw MetamodelException("NuTo::Metamodel::Build - number of outputs must be positive - set training data first.");
    
    if (mSupportPoints.GetDimInput()<=0)
        throw MetamodelException("NuTo::Metamodel::Build - number of inputs must be positive - set training data first.");

    if (!mSupportPoints.IsTransformationBuild())
        mSupportPoints.BuildTransformation();
        
    BuildDerived();
}

void Metamodel::Info()const
{
	mSupportPoints.Info();
}
