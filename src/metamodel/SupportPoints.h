// $Id$

/*******************************************************************************
Bauhaus-Universität Weimar
Author: Joerg F. Unger,  Septermber 2009
*******************************************************************************/


#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION


#include "math/FullMatrix_Def.h"
#include <boost/ptr_container/ptr_list.hpp>

namespace NuTo
{

class Transformation;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief stores the support points
class SupportPoints
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
	//! @brief constructor
    SupportPoints();
    //! @brief destructor
    ~SupportPoints();

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief clear support points
    void Clear();

    //! @brief info about support points
    void Info()const;

    //! @brief get number of support points
    inline int GetNumSupportPoints()const
	{
        return mSPOrigInput.cols();
	}

    //! @brief get input dimension of support points (input)
    inline int GetDimInput()const                              
	{
        return mSPOrigInput.rows();
	}

    //! @brief get output dimension of support points (input)
    inline int GetDimOutput()const
	{
        return mSPOrigOutput.rows();
	}

    //! @brief returns the input of the support points in a matrix
    inline const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetOrigSupportPointsInput()const
	{
	    return mSPOrigInput;
	}
    
    //! @brief returns the output of the support points in a matrix
    inline const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetOrigSupportPointsOutput()const
	{
	    return mSPOrigOutput;
	}
    
    //! @brief returns the input of the transformed support points in a matrix
    inline const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetTransformedSupportPointsInput()const
	{
	    return mSPTransInput;
	}
    
    //! @brief returns the input of the transformed support points in a matrix
    inline const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& GetTransformedSupportPointsOutput()const
	{
	    return mSPTransOutput;
	}
    
    //! @brief append a transformation for input points
    //! @brief in general, the order is orig_input -> outputtrans1 -> outputtrans2 -> trans_input
    void AppendTransformationInput(Transformation* rTransformation);

    //! @brief append a transformation for output points
    //! @brief in general, order is orig_output -> outputtrans1 -> outputtrans2 -> trans_output
    void AppendTransformationOutput(Transformation* rTransformation);

    //! @brief checks, if the transformation has been build
	bool IsTransformationBuild()const
	{
	    return mTransformationBuild;
	}

    //! @brief build transformation for support points
    void BuildTransformation();

    //! @brief set support points
    void SetSupportPoints(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSPOrigInput, const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rSPOrigOutput);

    //! @brief return the Weight vector as Matrix (use only from Python level, since everything is copied)
    inline const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> GetWeights(int rSample)const
    {
        return FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>(mWeight.size(),1,mWeight);
    }

    //! @brief Removes all applied Transformations for the input and the output points
    void ClearTransformations();
    
#ifndef SWIG
    inline double GetWeight(int rSample)const
	{
        return mWeight[rSample];
	}

    //! @brief perform forward transformation for inputs (from orig to transformed)
    //! @brief in general, the forward direction is orig_input -> inputtrans1 ->inputtrans2 -> metamodel -> outputtrans2 -> outputtrans1 -> orig_output
    void TransformForwardInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;
    
    //! @brief perform backward transformation for inputs (from transformed to orig
    void TransformBackwardInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;
    
    //! @brief perform forward transformation for outputs (from transformed to orig)
    //! @brief attention, this is exactly the backwards order, since transformations for outputs are given in revers order
    void TransformForwardOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;
    
    //! @brief perform backward transformation for outputs
    void TransformBackwardOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCoordinates)const;
    
#endif

    //! @brief calculate the mean value of the original inputs
    //! @param rMean ... vector of mean values (output)
    void GetMeanValueOriginalInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMean)const;

    //! @brief calculate the mean value of the transformed inputs
    //! @param rMean ... vector of mean values (output)
    void GetMeanValueTransformedInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMean)const;

    //! @brief calculate the mean value of the original outputs
    //! @param rMean ... vector of mean values (output)
    void GetMeanValueOriginalOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMean)const;

    //! @brief calculate the mean value of the transformed outputs
    //! @param rMean ... vector of mean values (output)
    void GetMeanValueTransformedOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMean)const;

    //! @brief calculate the variance of the original inputs
    //! @param rVariance ... vector of variances (output)
    void GetVarianceOriginalInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rVariance)const;

    //! @brief calculate the variance of the transformed inputs
    //! @param rVariance ... vector of variances (output)
    void GetVarianceTransformedInput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rVariance)const;

    //! @brief calculate the variance of the transformed outputs
    //! @param rVariance ... vector of variances (output)
    void GetVarianceOriginalOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rVariance)const;

    //! @brief calculate the variance of the original outputs
    //! @param rVariance ... vector of variances (output)
    void GetVarianceTransformedOutput(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rVariance)const;

    //! @brief calculate the covariance matrix of original inputs and outputs
    //! @param rCovarianceMatrix ... covariance matrix
    void GetCovarianceMatrixOriginal(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

    //! @brief calculate the covariance matrix of transformed inputs and outputs
    //! @param rCovarianceMatrix ... covariance matrix
    void GetCovarianceMatrixTransformed(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

    //! @brief calculate Pearson's correlation matrix
    //! @param rCorrelationMatrix  ... Pearson's correlation matrix using original support point coordinates
    void GetPearsonCorrelationMatrixOriginal(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix using transformed support point coordinates
    //! @param rCorrelationMatrix  ... Pearsons correlation matrix
    void GetPearsonCorrelationMatrixTransformed(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const;

    //! @brief calculate the confidence interval on Pearson's correlation coefficient using original support point coordinates
    //! @param rCorrelationMatrix  ... Pearson's correlation matrix using original support point coordinates
	//! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rMaxCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetPearsonCorrelationMatrixConfidenceIntervalsOriginal(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha = 0.05) const;

    //! @brief calculate the confidence interval on Pearson's correlation coefficient using transformed support point coordinates
    //! @param rCorrelationMatrix  ... Pearson's correlation matrix using original support point coordinates
	//! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rMaxCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetPearsonCorrelationMatrixConfidenceIntervalsTransformed(FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha = 0.05) const;
private:
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  mSPOrigInput;      //!< original inputs, each sample after another
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  mSPOrigOutput;     //!< original outputs, each sample after another
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  mSPTransInput;     //!< transformed inputs, each sample after another
    FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>  mSPTransOutput;    //!< transformed outputs, each sample after another

    std::vector<double> mWeight;    //!< weight of each support point
    
    boost::ptr_list<Transformation> mlTransformationInput;   //!< list of Transformations for inputs
    boost::ptr_list<Transformation> mlTransformationOutput;  //!< list of Transformations for outputs
	
	bool mTransformationBuild;
    
	//! @brief calculate the mean values of given data
	//! @param rData ... matrix of data points, each sample after another
	//! @param rMean ... vector of mean values (output)
	void CalculateMeanValues(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rData, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMean) const;

	//! @brief calculate the variance in data
	//! @param rData ... matrix of data points, each sample after another
	//! @param rVariance ... vector of variances
	void CalculateVariance(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rData, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rVariance) const;

	//! @brief calculate the covariance matrix
	//! @param rInputData ... matrix of input data points, each sample after another
	//! @param rOutputData ... matrix of output data points, each sample after another
	//! @param rCovarianceMatrix ... covariance matrix
	void CalculateCovarianceMatrix(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputData, const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputData, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

	//! @brief calculate Pearson's correlation matrix
	//! @param rInputData ... matrix of input data points, each sample after another
	//! @param rOutputData ... matrix of output data points, each sample after another
	//! @param rCorrelationMatrix ... Pearson's correlation matrix
	void CalculatePearsonCorrelationMatrix(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputData, const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputData, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const;

	//! @brief calculate Pearson's correlation matrix
	//! @param rInputData ... matrix of input data points, each sample after another
	//! @param rOutputData ... matrix of output data points, each sample after another
	//! @param rCorrelationMatrix ... Pearson's correlation matrix
	//! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rAlpha ... the confidence level is defined as (1-rAlpha)
	void CalculatePearsonCorrelationMatrixConfidenceIntervals(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputData, const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputData, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha) const;

};
} // namespace nuto
