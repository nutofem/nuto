#pragma once

#include <boost/ptr_container/ptr_list.hpp>
#include <Eigen/Core>
#include <vector>

namespace NuTo
{

class Transformation;

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief stores the support points
class SupportPoints
{

public:
    //! @brief constructor
    SupportPoints();
    //! @brief destructor
    ~SupportPoints();

    //! @brief clear support points
    void Clear();

    //! @brief info about support points
    void Info() const;

    //! @brief get number of support points
    inline int GetNumSupportPoints() const
    {
        return mSPOrigInput.cols();
    }

    //! @brief get input dimension of support points (input)
    inline int GetDimInput() const
    {
        return mSPOrigInput.rows();
    }

    //! @brief get output dimension of support points (input)
    inline int GetDimOutput() const
    {
        return mSPOrigOutput.rows();
    }

    //! @brief returns the input of the support points in a matrix
    inline const Eigen::MatrixXd& GetOrigSupportPointsInput() const
    {
        return mSPOrigInput;
    }

    //! @brief returns the output of the support points in a matrix
    inline const Eigen::MatrixXd& GetOrigSupportPointsOutput() const
    {
        return mSPOrigOutput;
    }

    //! @brief returns the input of the transformed support points in a matrix
    inline const Eigen::MatrixXd& GetTransformedSupportPointsInput() const
    {
        return mSPTransInput;
    }

    //! @brief returns the input of the transformed support points in a matrix
    inline const Eigen::MatrixXd& GetTransformedSupportPointsOutput() const
    {
        return mSPTransOutput;
    }

    //! @brief append a transformation for input points
    //! in general, the order is orig_input -> outputtrans1 -> outputtrans2 -> trans_input
    void AppendTransformationInput(Transformation* rTransformation);

    //! @brief append a transformation for output points
    //! in general, order is orig_output -> outputtrans1 -> outputtrans2 -> trans_output
    void AppendTransformationOutput(Transformation* rTransformation);

    //! @brief checks, if the transformation has been build
    bool IsTransformationBuild() const
    {
        return mTransformationBuild;
    }

    //! @brief build transformation for support points
    void BuildTransformation();

    //! @brief set support points
    void SetSupportPoints(const Eigen::MatrixXd& rSPOrigInput, const Eigen::MatrixXd& rSPOrigOutput);

    //! @brief return the Weight vector
    std::vector<double> GetWeights() const
    {
        return mWeight;
    }

    //! @brief Removes all applied Transformations for the input and the output points
    void ClearTransformations();

#ifndef SWIG
    inline double GetWeight(int rSample) const
    {
        return mWeight[rSample];
    }

    //! @brief perform forward transformation for inputs (from orig to transformed)
    //! in general, the forward direction is orig_input -> inputtrans1 ->inputtrans2 -> metamodel -> outputtrans2
    //! -> outputtrans1 -> orig_output
    void TransformForwardInput(Eigen::MatrixXd& rCoordinates) const;

    //! @brief perform backward transformation for inputs (from transformed to orig
    void TransformBackwardInput(Eigen::MatrixXd& rCoordinates) const;

    //! @brief perform forward transformation for outputs (from transformed to orig)
    //! @note attention, this is exactly the backwards order, since transformations for outputs are given in reverse
    //! order
    void TransformForwardOutput(Eigen::MatrixXd& rCoordinates) const;

    //! @brief perform backward transformation for outputs
    void TransformBackwardOutput(Eigen::MatrixXd& rCoordinates) const;

#endif

    //! @brief calculate the mean value of the original inputs
    //! @param rMean vector of mean values (output)
    void GetMeanValueOriginalInput(Eigen::MatrixXd& rMean) const;

    //! @brief calculate the mean value of the transformed inputs
    //! @param rMean vector of mean values (output)
    void GetMeanValueTransformedInput(Eigen::MatrixXd& rMean) const;

    //! @brief calculate the mean value of the original outputs
    //! @param rMean vector of mean values (output)
    void GetMeanValueOriginalOutput(Eigen::MatrixXd& rMean) const;

    //! @brief calculate the mean value of the transformed outputs
    //! @param rMean vector of mean values (output)
    void GetMeanValueTransformedOutput(Eigen::MatrixXd& rMean) const;

    //! @brief calculate the variance of the original inputs
    //! @param rVariance vector of variances (output)
    void GetVarianceOriginalInput(Eigen::MatrixXd& rVariance) const;

    //! @brief calculate the variance of the transformed inputs
    //! @param rVariance vector of variances (output)
    void GetVarianceTransformedInput(Eigen::MatrixXd& rVariance) const;

    //! @brief calculate the variance of the transformed outputs
    //! @param rVariance vector of variances (output)
    void GetVarianceOriginalOutput(Eigen::MatrixXd& rVariance) const;

    //! @brief calculate the variance of the original outputs
    //! @param rVariance vector of variances (output)
    void GetVarianceTransformedOutput(Eigen::MatrixXd& rVariance) const;

    //! @brief calculate the covariance matrix of original inputs and outputs
    //! @param rCovarianceMatrix covariance matrix
    void GetCovarianceMatrixOriginal(Eigen::MatrixXd& rCovarianceMatrix) const;

    //! @brief calculate the covariance matrix of transformed inputs and outputs
    //! @param rCovarianceMatrix covariance matrix
    void GetCovarianceMatrixTransformed(Eigen::MatrixXd& rCovarianceMatrix) const;

    //! @brief calculate Pearson's correlation matrix
    //! @param rCorrelationMatrix  Pearson's correlation matrix using original support point coordinates
    void GetPearsonCorrelationMatrixOriginal(Eigen::MatrixXd& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix using transformed support point coordinates
    //! @param rCorrelationMatrix  Pearsons correlation matrix
    void GetPearsonCorrelationMatrixTransformed(Eigen::MatrixXd& rCorrelationMatrix) const;

    //! @brief calculate the confidence interval on Pearson's correlation coefficient using original support point
    //! coordinates
    //! @param rCorrelationMatrix  Pearson's correlation matrix using original support point coordinates
    //! @param rMinCorrelationMatrix lower bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rMaxCorrelationMatrix upper bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rAlpha the confidence level is defined as (1-rAlpha)
    void GetPearsonCorrelationMatrixConfidenceIntervalsOriginal(Eigen::MatrixXd& rCorrelationMatrix,
                                                                Eigen::MatrixXd& rMinCorrelationMatrix,
                                                                Eigen::MatrixXd& rMaxCorrelationMatrix,
                                                                double rAlpha = 0.05) const;

    //! @brief calculate the confidence interval on Pearson's correlation coefficient using transformed support point
    //! coordinates
    //! @param rCorrelationMatrix  Pearson's correlation matrix using original support point coordinates
    //! @param rMinCorrelationMatrix lower bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rMaxCorrelationMatrix upper bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rAlpha the confidence level is defined as (1-rAlpha)
    void GetPearsonCorrelationMatrixConfidenceIntervalsTransformed(Eigen::MatrixXd& rCorrelationMatrix,
                                                                   Eigen::MatrixXd& rMinCorrelationMatrix,
                                                                   Eigen::MatrixXd& rMaxCorrelationMatrix,
                                                                   double rAlpha = 0.05) const;

private:
    Eigen::MatrixXd mSPOrigInput; //!< original inputs, each sample after another
    Eigen::MatrixXd mSPOrigOutput; //!< original outputs, each sample after another
    Eigen::MatrixXd mSPTransInput; //!< transformed inputs, each sample after another
    Eigen::MatrixXd mSPTransOutput; //!< transformed outputs, each sample after another

    std::vector<double> mWeight; //!< weight of each support point

    boost::ptr_list<Transformation> mlTransformationInput; //!< list of Transformations for inputs
    boost::ptr_list<Transformation> mlTransformationOutput; //!< list of Transformations for outputs

    bool mTransformationBuild;

    //! @brief calculate the mean values of given data
    //! @param rData matrix of data points, each sample after another
    //! @param rMean vector of mean values (output)
    void CalculateMeanValues(const Eigen::MatrixXd& rData, Eigen::MatrixXd& rMean) const;

    //! @brief calculate the variance in data
    //! @param rData matrix of data points, each sample after another
    //! @param rVariance vector of variances
    void CalculateVariance(const Eigen::MatrixXd& rData, Eigen::MatrixXd& rVariance) const;

    //! @brief calculate the covariance matrix
    //! @param rInputData matrix of input data points, each sample after another
    //! @param rOutputData matrix of output data points, each sample after another
    //! @param rCovarianceMatrix covariance matrix
    void CalculateCovarianceMatrix(const Eigen::MatrixXd& rInputData, const Eigen::MatrixXd& rOutputData,
                                   Eigen::MatrixXd& rCovarianceMatrix) const;

    //! @brief calculate Pearson's correlation matrix
    //! @param rInputData matrix of input data points, each sample after another
    //! @param rOutputData matrix of output data points, each sample after another
    //! @param rCorrelationMatrix Pearson's correlation matrix
    void CalculatePearsonCorrelationMatrix(const Eigen::MatrixXd& rInputData, const Eigen::MatrixXd& rOutputData,
                                           Eigen::MatrixXd& rCorrelationMatrix) const;

    //! @brief calculate Pearson's correlation matrix
    //! @param rInputData matrix of input data points, each sample after another
    //! @param rOutputData matrix of output data points, each sample after another
    //! @param rCorrelationMatrix Pearson's correlation matrix
    //! @param rMinCorrelationMatrix lower bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rMaxCorrelationMatrix upper bound of the confidence interval on the coefficients of Pearson's
    //! correlation matrix
    //! @param rAlpha the confidence level is defined as (1-rAlpha)
    void CalculatePearsonCorrelationMatrixConfidenceIntervals(
            const Eigen::MatrixXd& rInputData, const Eigen::MatrixXd& rOutputData, Eigen::MatrixXd& rCorrelationMatrix,
            Eigen::MatrixXd& rMinCorrelationMatrix, Eigen::MatrixXd& rMaxCorrelationMatrix, double rAlpha) const;
};
} // namespace nuto
