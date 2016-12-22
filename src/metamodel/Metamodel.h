// $Id$

#pragma once
#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif  // ENABLE_SERIALIZATION


#include <random>
#include <vector>

#include "metamodel/MetamodelException.h"
#include "base/NuToObject.h"
#include "metamodel/SupportPoints.h"


namespace NuTo
{
template<class T, int rows, int cols> class FullMatrix;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all metamodels in NuTo
class Metamodel : public virtual NuToObject
{
#ifdef ENABLE_SERIALIZATION
	friend class boost::serialization::access;
#endif // ENABLE_SERIALIZATION

public:
	//! @brief constructor
    Metamodel();
    
#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    void Build();
    virtual void BuildDerived()=0;
    
    void AppendMinMaxTransformationInput(double min, double max);
    void AppendMinMaxTransformationInput(int coordinate, double min, double max);
    void AppendMinMaxTransformationOutput(double min, double max);
    void AppendMinMaxTransformationOutput(int coordinate, double min, double max);

    //! @brief ... append for all coordinates of the input points a transformation to zero mean and unit variance to the list of transformations.
    void AppendZeroMeanUnitVarianceTransformationInput();
    
    //! @brief ... append for a specific coordinate of the input points a transformation to zero mean and unit variance to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate);

    //! @brief ... append for all coordinates of the output points a transformation to zero mean and unit variance to the list of transformations.
    void AppendZeroMeanUnitVarianceTransformationOutput();
    
    //! @brief ... append for a specific coordinate of the output points a transformation to zero mean and unit variance to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate);

    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetOriginalSupportPointsInput()const;
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetOriginalSupportPointsOutput()const;
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetTransformedSupportPointsInput()const;
    FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> GetTransformedSupportPointsOutput()const;
    void SetSupportPoints(int rDimInput, int rDimOutput, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rInputCoordinates, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> rOutputCoordinates);
    void BuildTransformation();
    void InitRandomNumberGenerator(int rSeed);
    double RandomDouble();
    virtual void Info()const;
    void Solve(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const;
    void SolveConfidenceInterval(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates,
                                                    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMax)const;
    virtual void SolveTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const=0;
    virtual void SolveConfidenceIntervalTransformed(const FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates,
                                            NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinatesMax)const
    {
        throw MetamodelException("Metamodel::SolveConfidenceIntervalTransformed - not implemented for this kind of metamodel.");
    }

    //! @brief calculate the sample mean for each support point using original support point coordinates
    //! @param rInputMean ... vector of mean values of the input support point coordinates
    //! @param rOutputMean ... vector of mean values of the output support point coordinates
    void GetOriginalSupportPointsMeanValue(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean) const;

    //! @brief calculate the sample mean for each support point using transformed support point coordinates
    //! @param rInputMean ... vector of mean values of the input support point coordinates
    //! @param rOutputMean ... vector of mean values of the output support point coordinates
    void GetTransformedSupportPointsMeanValue(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean) const;

    //! @brief calculate the sample variance for each support point using original support point coordinates
    //! @param rInputVariance ... vector of sample variances of the input support point coordinates
    //! @param rOutputVariance ... vector of sample variances of the output support point coordinates
    void GetOriginalSupportPointsVariance(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputVariance, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputVariance) const;

    //! @brief calculate the sample variance for each support point using transformed support point coordinates
    //! @param rInputVariance ... vector of sample variances of the input support point coordinates
    //! @param rOutputVariance ... vector of sample variances of the output support point coordinates
    void GetTransformedSupportPointsVariance(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputVariance, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputVariance) const;

    //! @brief calculate the covariance matrix of the support points using original support point coordinates
    //! @param rCovarianceMatrix ... covariance matrix
    void GetOriginalSupportPointsCovarianceMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

    //! @brief calculate the covariance matrix of the support points using transformed support point coordinates
    //! @param rCovarianceMatrix ... covariance matrix
    void GetTransformedSupportPointsCovarianceMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    void GetOriginalSupportPointsPearsonCorrelationMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    void GetTransformedSupportPointsPearsonCorrelationMatrix(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
	//! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetOriginalSupportPointsPearsonCorrelationMatrixConfidenceInterval(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha = 0.05) const;

    //! @brief calculate Pearson's correlation matrix of the support points using transformed support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
	//! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's correlation matrix
	//! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetTransformedSupportPointsPearsonCorrelationMatrixConfidenceInterval(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCorrelationMatrix, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rMinCorrelationMatrix, FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rMaxCorrelationMatrix, double rAlpha = 0.05) const;
protected:

    NuTo::SupportPoints mSupportPoints;
    std::mt19937_64 mRandomNumberGenerator;
};
} //namespace NuTo
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::Metamodel)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::Metamodel)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
