// $Id$

/*******************************************************************************
 Bauhaus-University Weimar
 Author: Joerg F. Unger ,  September 2009
*******************************************************************************/

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/ptr_container/serialize_ptr_list.hpp>
#endif //ENABLE_SERIALIZATION

#include "nuto/metamodel/SupportPoints.h"
#include "nuto/metamodel/MetamodelException.h"

NuTo::SupportPoints::SupportPoints()
{
    mTransformationBuild = false;
}

NuTo::SupportPoints::~SupportPoints()
{
    ClearTransformations();
}

//! @brief clear support points
void NuTo::SupportPoints::Clear()
{
    mSPOrigInput.Resize(0,0);
    mSPTransInput.Resize(0,0);
    mSPTransOutput.Resize(0,0);
    mWeight.resize(0);

    mlTransformationInput.clear();
    mlTransformationOutput.clear();

    mTransformationBuild = false;
}

//! @brief info about support points
void NuTo::SupportPoints::Info()const
{
    throw MetamodelException("[SupportPoints::Info()] not yet implemented.");
}

void NuTo::SupportPoints::BuildTransformation()
{
	mSPTransInput = mSPOrigInput;
	for (boost::ptr_list<Transformation>::iterator it = mlTransformationInput.begin(); it!=mlTransformationInput.end();it++)
	{
		it->Build(mSPTransInput);
        it->TransformForward(mSPTransInput);
    }
    
	mSPTransOutput = mSPOrigOutput;
	for (boost::ptr_list<Transformation>::iterator it = mlTransformationOutput.begin(); it!=mlTransformationOutput.end();it++)
	{
		it->Build(mSPTransOutput);
        it->TransformForward(mSPTransOutput);
    }
	mTransformationBuild = true;
}


void NuTo::SupportPoints::AppendTransformationInput(Transformation* rTransformation)
{
    mlTransformationInput.push_back(rTransformation);
}

void NuTo::SupportPoints::AppendTransformationOutput(Transformation* rTransformation)
{
    mlTransformationOutput.push_back(rTransformation);
}

void NuTo::SupportPoints::SetSupportPoints(const FullMatrix<double>& rSPOrigInput, const FullMatrix<double>& rSPOrigOutput)
{
    if (rSPOrigInput.GetNumColumns()!=rSPOrigOutput.GetNumColumns())
    {
        throw MetamodelException("[NuTo::SupportPoints::SetSupportPoints] Number of columns for input and output must be identical (=number of samples).");   
    }
    
    mSPOrigInput = rSPOrigInput;
    mSPOrigOutput  = rSPOrigOutput;

	mWeight.resize(GetNumSupportPoints(),1);
    
    // clear all transformations, since it is no longer for sure that the dimensions remain identical
    mTransformationBuild = false;
    ClearTransformations();
}

//! @brief perform forward transformation for inputs (from orig to transformed)
void NuTo::SupportPoints::TransformForwardInput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_iterator it=mlTransformationInput.begin();it!=mlTransformationInput.end(); it++)
        it->TransformForward(rCoordinates);
}

//! @brief perform backward transformation for inputs  (from transformed to orig)
void NuTo::SupportPoints::TransformBackwardInput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_reverse_iterator it=mlTransformationInput.rbegin();it!=mlTransformationInput.rend(); it++)
        it->TransformBackward(rCoordinates);
}

//! @brief perform forward transformation for outputs (from transformed to orig)
//! @brief attention, this is exactly the backwards order, since transformations for outputs are given in revers order
//! @brief orig->trans1->trans2-transformed
void NuTo::SupportPoints::TransformForwardOutput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_reverse_iterator it=mlTransformationOutput.rbegin();it!=mlTransformationOutput.rend(); it++)
        it->TransformBackward(rCoordinates);
}    

//! @brief perform backward transformation for outputs
//! @brief attention, this is exactly the backwards order, since transformations for given in revers order
void NuTo::SupportPoints::TransformBackwardOutput(FullMatrix<double>& rCoordinates)const
{
    for (boost::ptr_list<Transformation>::const_iterator it=mlTransformationOutput.begin();it!=mlTransformationOutput.end(); it++)
        it->TransformForward(rCoordinates);
}

//! @brief Clears all the transformations for input and output
void NuTo::SupportPoints::ClearTransformations()
{
    mlTransformationInput.clear();
    mlTransformationOutput.clear();
}

// calculate mean values of original inputs
void NuTo::SupportPoints::GetMeanValueOriginalInput(FullMatrix<double>& rMean)const
{
	this->CalculateMeanValues(this->mSPOrigInput, rMean);
}

// calculate mean values of original outputs
void NuTo::SupportPoints::GetMeanValueOriginalOutput(FullMatrix<double>& rMean)const
{
	this->CalculateMeanValues(this->mSPOrigOutput, rMean);
}

// calculate variance of original inputs
void NuTo::SupportPoints::GetVarianceOriginalInput(FullMatrix<double>& rVariance)const
{
	this->CalculateVariance(this->mSPOrigInput, rVariance);
}

// calculate mean values of original outputs
void NuTo::SupportPoints::GetVarianceOriginalOutput(FullMatrix<double>& rVariance)const
{
	this->CalculateVariance(this->mSPOrigOutput, rVariance);
}

// calculate covariance matrix
void NuTo::SupportPoints::GetCovarianceMatrixOriginal(FullMatrix<double>& rCovarianceMatrix) const
{
	this->CalculateCovarianceMatrix(this->mSPOrigInput, this->mSPOrigOutput, rCovarianceMatrix);
}

// calculate Pearson's correlation matrix
void NuTo::SupportPoints::GetPearsonCorrelationMatrixOriginal(FullMatrix<double>& rCorrelationMatrix) const
{
	this->CalculatePearsonCorrelationMatrix(this->mSPOrigInput, this->mSPOrigOutput, rCorrelationMatrix);
}

// calculate mean values
void NuTo::SupportPoints::CalculateMeanValues(const FullMatrix<double>& rData, FullMatrix<double>& rMean) const
{
	int numRows = rData.GetNumRows();
	int numSamples = rData.GetNumColumns();
	if(numSamples < 1)
	{
		throw MetamodelException("[NuTo::SupportPoints::CalculateMeanValues] number of samples must be larger than zero.");
	}
	double factor = 1.0/static_cast<double>(numSamples);
	rMean.Resize(numRows,1);

	// calculate mean value
	for(int rowCount = 0; rowCount < numRows; rowCount++)
	{
		double mean=0.0;
	    const double *dataPtr = &rData.mEigenMatrix.data()[rowCount];
	    for (int sample=0; sample<numSamples; sample++)
		{
	        mean += factor * (*dataPtr);
	        dataPtr+=numRows;
		}
	    rMean.SetValue(rowCount,0,mean);
	}
}

// calculate variances
void NuTo::SupportPoints::CalculateVariance(const FullMatrix<double>& rData, FullMatrix<double>& rVariance) const
{
	int numRows = rData.GetNumRows();
	int numSamples = rData.GetNumColumns();
	if(numSamples < 2)
	{
		throw MetamodelException("[NuTo::SupportPoints::CalculateVariance] number of samples must be larger than one.");
	}
	rVariance.Resize(numRows,1);
	double factor = 1.0/static_cast<double>(numSamples-1);

	// calculate mean values
	FullMatrix<double> meanVector;
	this->CalculateMeanValues(rData, meanVector);

	// calculate variance
	for(int rowCount = 0; rowCount < numRows; rowCount++)
	{
		double variance=0.0;
		double mean=meanVector.GetValue(rowCount,0);
	    const double *dataPtr = &rData.mEigenMatrix.data()[rowCount];
	    for (int sample=0; sample<numSamples; sample++)
		{
	    	double delta = *dataPtr - mean;
	        variance += factor * delta * delta;
	        dataPtr+=numRows;
		}
	    rVariance.SetValue(rowCount,0,variance);
	}
}

// calculate covariance matrix
void NuTo::SupportPoints::CalculateCovarianceMatrix(const FullMatrix<double>& rInputData, const FullMatrix<double>& rOutputData, FullMatrix<double>& rCovarianceMatrix) const
{
	// get data
	int numInputData = rInputData.GetNumRows();
	int numOutputData = rOutputData.GetNumRows();
	int numSamples = rInputData.GetNumColumns();
	if(numSamples != rOutputData.GetNumColumns())
	{
		throw MetamodelException("[NuTo::SupportPoints::CalculateCovarianceMatrix] number of samples in input data and number of samples in output data must be equal.");
	}
	if(numSamples < 2)
	{
		throw MetamodelException("[NuTo::SupportPoints::CalculateCovarianceMatrix] number of samples must be larger than one.");
	}

	// calculate mean values
	FullMatrix<double> inputDataMeanVector;
	this->CalculateMeanValues(rInputData, inputDataMeanVector);
	FullMatrix<double> outputDataMeanVector;
	this->CalculateMeanValues(rOutputData, outputDataMeanVector);

	// calculate covariance matrix
	rCovarianceMatrix.Resize(numInputData+numOutputData,numInputData+numOutputData);
	for(int row = 0; row < numInputData; row++)
	{
		double meanRow = inputDataMeanVector.GetValue(row,0);
		// covariance input - input
		for(int column = row; column < numInputData; column++)
		{
			double cov = 0;
			double meanColumn = inputDataMeanVector.GetValue(column,0);
			const double *dataRowPtr = &rInputData.mEigenMatrix.data()[row];
			const double *dataColumnPtr = &rInputData.mEigenMatrix.data()[column];
			for(int sample = 0; sample < numSamples; sample++)
			{
				cov += (*dataRowPtr - meanRow) * (*dataColumnPtr - meanColumn);
				dataRowPtr += numInputData;
				dataColumnPtr += numInputData;
			}
			cov /= static_cast<double>(numSamples-1);
			rCovarianceMatrix.SetValue(row,column, cov);
			if(row != column)
			{
				rCovarianceMatrix.SetValue(column, row, cov);
			}
		}
		// covariance input - output / output - input
		for(int column = numInputData; column < numInputData + numOutputData; column++)
		{
			double cov = 0;
			double meanColumn = outputDataMeanVector.GetValue(column - numInputData,0);
			const double *dataRowPtr = &rInputData.mEigenMatrix.data()[row];
			const double *dataColumnPtr = &rOutputData.mEigenMatrix.data()[column - numInputData];
			for(int sample = 0; sample < numSamples; sample++)
			{
				cov += (*dataRowPtr - meanRow) * (*dataColumnPtr - meanColumn);
				dataRowPtr += numInputData;
				dataColumnPtr += numOutputData;
			}
			cov /= static_cast<double>(numSamples-1);
			rCovarianceMatrix.SetValue(row,column, cov);
			rCovarianceMatrix.SetValue(column, row, cov);
		}
	}
	// covariance output - output
	for(int row = numInputData; row < numInputData + numOutputData; row++)
	{
		double meanRow = outputDataMeanVector.GetValue(row - numInputData,0);
		for(int column = row; column < numInputData + numOutputData; column++)
		{
			double cov = 0;
			double meanColumn = outputDataMeanVector.GetValue(column - numInputData,0);
			const double *dataRowPtr = &rOutputData.mEigenMatrix.data()[row - numInputData];
			const double *dataColumnPtr = &rOutputData.mEigenMatrix.data()[column - numInputData];
			for(int sample = 0; sample < numSamples; sample++)
			{
				cov += (*dataRowPtr - meanRow) * (*dataColumnPtr - meanColumn);
				dataRowPtr += numOutputData;
				dataColumnPtr += numOutputData;
			}
			cov /= static_cast<double>(numSamples-1);
			rCovarianceMatrix.SetValue(row,column, cov);
			if(row != column)
			{
				rCovarianceMatrix.SetValue(column, row, cov);
			}
		}
	}
}

void NuTo::SupportPoints::CalculatePearsonCorrelationMatrix(const FullMatrix<double>& rInputData, const FullMatrix<double>& rOutputData, FullMatrix<double>& rCorrelationMatrix) const
{
	// calculate covariance matrix
	this->CalculateCovarianceMatrix(rInputData, rOutputData, rCorrelationMatrix);

	// calculate standard deviation of diagonal values
	int dim = rCorrelationMatrix.GetNumRows();
	std::vector<double> stddev(dim);
	assert(dim == rCorrelationMatrix.GetNumColumns());
	assert(dim == rInputData.GetNumRows() + rOutputData.GetNumRows());
	for(int dimCount = 0; dimCount < dim; dimCount++)
	{
		assert(rCorrelationMatrix.GetValue(dimCount,dimCount) > 0.0);
		stddev[dimCount] = sqrt(rCorrelationMatrix.GetValue(dimCount,dimCount));
		std::cout << "stddev: " << stddev[dimCount] << std::endl;
	}

	// scale covariance matrix
	for(int row = 0; row < dim; row++)
	{
		for(int column = 0; column < dim; column++)
		{
			rCorrelationMatrix.SetValue(row,column, rCorrelationMatrix.GetValue(row,column)/(stddev[row]*stddev[column]));
		}
	}
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::SupportPoints::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::SupportPoints::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::SupportPoints::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::SupportPoints::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::SupportPoints::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::SupportPoints::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::SupportPoints::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize SupportPoints" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(mSPOrigInput)
       & BOOST_SERIALIZATION_NVP(mSPOrigOutput)
       & BOOST_SERIALIZATION_NVP(mSPTransInput)
       & BOOST_SERIALIZATION_NVP(mSPTransOutput)
       & BOOST_SERIALIZATION_NVP(mWeight)
       & BOOST_SERIALIZATION_NVP(mlTransformationInput)
       & BOOST_SERIALIZATION_NVP(mlTransformationOutput)
       & BOOST_SERIALIZATION_NVP(mTransformationBuild);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize SupportPoints" << std::endl;
#endif
}
#endif // ENABLE_SERIALIZATION
