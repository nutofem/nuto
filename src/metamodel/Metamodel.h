#pragma once


#include <random>
#include <vector>

#include "base/Exception.h"
#include "metamodel/SupportPoints.h"


namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief ... standard abstract class for all metamodels in NuTo
class Metamodel
{

public:
    //! @brief constructor
    Metamodel();
    virtual ~Metamodel() = default;

    void Build();
    virtual void BuildDerived() = 0;

    void AppendMinMaxTransformationInput(double min, double max);
    void AppendMinMaxTransformationInput(int coordinate, double min, double max);
    void AppendMinMaxTransformationOutput(double min, double max);
    void AppendMinMaxTransformationOutput(int coordinate, double min, double max);

    //! @brief ... append for all coordinates of the input points a transformation to zero mean and unit variance to the
    //! list of transformations.
    void AppendZeroMeanUnitVarianceTransformationInput();

    //! @brief ... append for a specific coordinate of the input points a transformation to zero mean and unit variance
    //! to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationInput(int rCoordinate);

    //! @brief ... append for all coordinates of the output points a transformation to zero mean and unit variance to
    //! the list of transformations.
    void AppendZeroMeanUnitVarianceTransformationOutput();

    //! @brief ... append for a specific coordinate of the output points a transformation to zero mean and unit variance
    //! to the list of transformations.
    //! @param rCoordinate ... coordinate for which the transformation is performed
    void AppendZeroMeanUnitVarianceTransformationOutput(int rCoordinate);

    Eigen::MatrixXd GetOriginalSupportPointsInput() const;
    Eigen::MatrixXd GetOriginalSupportPointsOutput() const;
    Eigen::MatrixXd GetTransformedSupportPointsInput() const;
    Eigen::MatrixXd GetTransformedSupportPointsOutput() const;
    void SetSupportPoints(int rDimInput, int rDimOutput, Eigen::MatrixXd rInputCoordinates,
                          Eigen::MatrixXd rOutputCoordinates);
    void BuildTransformation();
    void InitRandomNumberGenerator(int rSeed);
    double RandomDouble();
    virtual void Info() const;
    void Solve(const Eigen::MatrixXd& rInputCoordinates, Eigen::MatrixXd& rOutputCoordinates) const;
    void SolveConfidenceInterval(const Eigen::MatrixXd& rInputCoordinates, Eigen::MatrixXd& rOutputCoordinates,
                                 Eigen::MatrixXd& rOutputCoordinatesMin, Eigen::MatrixXd& rOutputCoordinatesMax) const;
    virtual void SolveTransformed(const Eigen::MatrixXd& rInputCoordinates,
                                  Eigen::MatrixXd& rOutputCoordinates) const = 0;
    virtual void SolveConfidenceIntervalTransformed(const Eigen::MatrixXd& rInputCoordinates,
                                                    Eigen::MatrixXd& rOutputCoordinates,
                                                    Eigen::MatrixXd& rOutputCoordinatesMin,
                                                    Eigen::MatrixXd& rOutputCoordinatesMax) const
    {
        throw Exception("Metamodel::SolveConfidenceIntervalTransformed - not implemented for this kind of metamodel.");
    }

    //! @brief calculate the sample mean for each support point using original support point coordinates
    //! @param rInputMean ... vector of mean values of the input support point coordinates
    //! @param rOutputMean ... vector of mean values of the output support point coordinates
    void GetOriginalSupportPointsMeanValue(Eigen::MatrixXd& rInputMean, Eigen::MatrixXd& rOutputMean) const;

    //! @brief calculate the sample mean for each support point using transformed support point coordinates
    //! @param rInputMean ... vector of mean values of the input support point coordinates
    //! @param rOutputMean ... vector of mean values of the output support point coordinates
    void GetTransformedSupportPointsMeanValue(Eigen::MatrixXd& rInputMean, Eigen::MatrixXd& rOutputMean) const;

    //! @brief calculate the sample variance for each support point using original support point coordinates
    //! @param rInputVariance ... vector of sample variances of the input support point coordinates
    //! @param rOutputVariance ... vector of sample variances of the output support point coordinates
    void GetOriginalSupportPointsVariance(Eigen::MatrixXd& rInputVariance, Eigen::MatrixXd& rOutputVariance) const;

    //! @brief calculate the sample variance for each support point using transformed support point coordinates
    //! @param rInputVariance ... vector of sample variances of the input support point coordinates
    //! @param rOutputVariance ... vector of sample variances of the output support point coordinates
    void GetTransformedSupportPointsVariance(Eigen::MatrixXd& rInputVariance, Eigen::MatrixXd& rOutputVariance) const;

    //! @brief calculate the covariance matrix of the support points using original support point coordinates
    //! @param rCovarianceMatrix ... covariance matrix
    void GetOriginalSupportPointsCovarianceMatrix(Eigen::MatrixXd& rCovarianceMatrix) const;

    //! @brief calculate the covariance matrix of the support points using transformed support point coordinates
    //! @param rCovarianceMatrix ... covariance matrix
    void GetTransformedSupportPointsCovarianceMatrix(Eigen::MatrixXd& rCovarianceMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    void GetOriginalSupportPointsPearsonCorrelationMatrix(Eigen::MatrixXd& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    void GetTransformedSupportPointsPearsonCorrelationMatrix(Eigen::MatrixXd& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix of the support points using original support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    //! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetOriginalSupportPointsPearsonCorrelationMatrixConfidenceInterval(Eigen::MatrixXd& rCorrelationMatrix,
                                                                            Eigen::MatrixXd& rMinCorrelationMatrix,
                                                                            Eigen::MatrixXd& rMaxCorrelationMatrix,
                                                                            double rAlpha = 0.05) const;

    //! @brief calculate Pearson's correlation matrix of the support points using transformed support point coordinates
    //! @param rCorrelationMatrix ... Pearson's correlation matrix
    //! @param rMinCorrelationMatrix ... lower bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rCorrelationMatrix ... upper bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rAlpha ... the confidence level is defined as (1-rAlpha)
    void GetTransformedSupportPointsPearsonCorrelationMatrixConfidenceInterval(Eigen::MatrixXd& rCorrelationMatrix,
                                                                               Eigen::MatrixXd& rMinCorrelationMatrix,
                                                                               Eigen::MatrixXd& rMaxCorrelationMatrix,
                                                                               double rAlpha = 0.05) const;

    void SetVerboseLevel(unsigned short verboseLevel);

protected:
    NuTo::SupportPoints mSupportPoints;
    std::mt19937_64 mRandomNumberGenerator;
    unsigned short mVerboseLevel;
};
} // namespace NuTo
