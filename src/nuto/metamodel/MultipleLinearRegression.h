// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

//parent
#include "nuto/metamodel/Metamodel.h"

namespace NuTo
{
//! @brief Class for the multiple linear regression model.
//! @author Stefan Eckardt
//! @date November 2010
/*!
 *  In the multiple linear regression model the relationship between the dependent variable \f$y\f$ and the \f$p\f$ regressor variables
 *  \f$\boldsymbol{x}^T = \left[ x_1, x_2, \ldots, x_p\right]\f$ is given by
 *  \f[ y = \beta_0 + \sum_{i=1}^p \beta_i x_i,\f]
 *  where \f$\beta_0\f$ and \f$\beta_i\f$ are the regression coefficients.
 *
 *  The implementation of this model is based on:
 *  - William W. Hines, Douglas C. Montgomery, David M. Goldsman and Connie M. Borror. Probability and Statistics in Engineering. 4-th edition. John Wiley & Sons, Inc. 2003.
 *  - Christian Bucher. Computational Analysis of Randomness in Structural Mechanics. CRC Press/Balkema. 2009.
 *
 */
class MultipleLinearRegression : public virtual NuTo::Metamodel
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION
public:
    //! @brief default constructor
    MultipleLinearRegression();

#ifdef ENABLE_SERIALIZATION
    //! @brief ... restore the object from a file
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
    //! @brief ... save the object to a file
    void Restore (const std::string &filename, std::string rType );

	//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
    //! @param filename ... filename
    //! @param aType ... type of file, either BINARY, XML or TEXT
	void Save (const std::string &filename, std::string rType )const;

	//! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

    //! @brief determine regression parameters
    void BuildDerived();

    //! @brief ... calculate approximation (in transformed space)
    //! @param rInputCoordinates ... matrix of input data points (transformed)
    //! @param rOutputCoordinates ... vector of output data (transformed)
    void SolveTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInputCoordinates, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputCoordinates)const;

    //! @brief ... calculate the residual between the support point data and the linear regression
    //! @param rResidualVector ... residual vector
    void GetSupportPointsResidual(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResidualVector) const;

    //! @brief ... calculate the residual between the transformed support point data and the linear regression
    //! @param rResidualVector ... residual vector (transformed)
    void GetSupportPointsResidualTransformed(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rResidualVector) const;

    //! @brief ... calculate the mean square error using original support points
    //! @return mean square error
    double GetMeanSquareError() const;

    //! @brief ... calculate the mean square error using transformed support points
    //! @return mean square error
    double GetMeanSquareErrorTransformed() const;

    //! @brief ... get regression coefficients
    //! @param rRegressionCoefficients ... vector of regression coefficients
    void GetRegressionCoefficientsTransformed(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRegressionCoefficients) const;

    //! @brief ... get covariance matrix of regression coefficients
    //! @param rCovarianceMatrix ... covariance matrix
    void GetRegressionCoefficientsCovarianceMatrixTransformed(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rCovarianceMatrix) const;

    //! @brief ... get the (1-rAlpha) confidence interval on the regression coefficients
    //! @param rRegressionCoefficients ... vector of regression coefficients
    //! @param rRegressionCoefficientsMin ... lower bound of the confidence interval on regression coefficients
    //! @param rRegressionCoefficientsMax ... upper bound of the confidence interval on regression coefficients
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void GetRegressionCoefficientsConfidenceIntervalsTransformed(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRegressionCoefficients, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRegressionCoefficientsMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rRegressionCoefficientsMax, double rAlpha = 0.05) const;

    //! @brief ... calculate the (1-rAlpha) confidence interval on the mean response of the linear regression
    //! @param rInput ... input point coordinates
    //! @param rOutputMean ... output of the regression model
    //! @param rOutputMin ... lower bound of the confidence interval
    //! @param rOutputMax ... upper bound of the confidence interval
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void SolveConfidenceIntervalsMeanResponse(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInput, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMax, double rAlpha = 0.05) const;

    //! @brief ... calculate the (1-rAlpha) confidence interval on predictions of the linear regression
    //! @param rInput ... input point coordinates
    //! @param rOutputMean ... output of the regression model
    //! @param rOutputMin ... lower bound of the confidence interval
    //! @param rOutputMax ... upper bound of the confidence interval
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void SolveConfidenceIntervalsPrediction(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInput, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMax, double rAlpha = 0.05) const;

    //! @brief calculate the coefficient of determination (correlation between the outputs predicted by the regression model and the actual data) using original support point coordinates
    //! @return ... coefficient of determination
    double GetCoefficientOfDetermination() const;

    //! @brief calculate the coefficient of determination (correlation between the outputs predicted by the regression model and the actual data) using transformed support point coordinates
    //! @return ... coefficient of determination
    double GetCoefficientOfDeterminationTransformed() const;

    //! @brief calculate the adjusted coefficient of determination (taking into account the sample size) using original support point data
    //! @return ... adjusted coefficient of determination
    double GetAdjustedCoefficientOfDetermination() const;

    //! @brief calculate the adjusted coefficient of determination (taking into account the sample size) using transformed support point data
    //! @return ... adjusted coefficient of determination
    double GetAdjustedCoefficientOfDeterminationTransformed() const;

    //! @brief ... get total sum of squares (original support point data)
    //! @return total sum of squares (original support point data)
    double GetTotalSumOfSquares() const;

    //! @brief ... get total sum of squares (transformed support point data)
    //! @return total sum of squares (transformed support point data)
    double GetTotalSumOfSquaresTransformed() const;

    //! @brief ... get error sum of squares (original support point data)
    //! @return error sum of squares (original support point data)
    double GetErrorSumOfSquares() const;

    //! @brief ... get error sum of squares (transformed support point data)
    //! @return error sum of squares (transformed support point data)
    double GetErrorSumOfSquaresTransformed() const;

    //! @brief ... get regression sum of squares (original support point data)
    //! @return regression sum of squares
    double GetRegressionSumOfSquares() const;

    //! @brief ... get regression sum of squares (transformed support point data)
    //! @return regression sum of squares
    double GetRegressionSumOfSquaresTransformed() const;

    //! @brief ... test for the significance of the regression
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    //! @return ... true if the regression is significant, false otherwise
    bool TestRegressionSignificance(double rAlpha=0.05) const;

    //! @brief ... test for the significance of the regression (transformed support points)
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    //! @return ... true if the regression is significant, false otherwise
    bool TestRegressionSignificanceTransformed(double rAlpha=0.05) const;

    //! @brief ... test for the significance of the individual regression coefficients
    //! @param rTestResult ... result of the significance test
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void TestRegressionCoefficientsSignificanceTransformed(std::vector<bool>& rTestResult, double rAlpha=0.05) const;

    //! @brief test the significance of a set of regressor variables/coefficients using original support point coordinates
    //! @param rTestCoefficients ... matrix of identifiers (number of identifiers, 1) defining the subset of regressor variables/coefficients (the intercept \f$\beta_0\f$ is not tested) for which the test is performed
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    //! @return true if the contribution of the subset is significant, false otherwise
    bool TestGeneralRegressionSignificance(const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& rTestCoefficients, double rAlpha=0.05) const;

    //! @brief test the significance of a set of regressor variables/coefficients using transformed support point coordinates
    //! @param rTestCoefficients ... matrix of identifiers (number of identifiers, 1) defining the subset of regressor variables/coefficients (the intercept \f$\beta_0\f$ is not tested) for which the test is performed
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    //! @return true if the contribution of the subset is significant, false otherwise
    bool TestGeneralRegressionSignificanceTransformed(const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& rTestCoefficients, double rAlpha=0.05) const;

	//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
    //!            in case of restoring from a file with the wrong object type, the file id is printed
    //! @return    class name
    std::string GetTypeId()const;
protected:
    //! @brief ... calculate the (1-rAlpha) confidence interval on the mean response of the linear regression
    //! @param rInput ... input point coordinates
    //! @param rOutputMean ... output of the regression model
    //! @param rOutputMin ... lower bound of the confidence interval
    //! @param rOutputMax ... upper bound of the confidence interval
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void SolveConfidenceIntervalsMeanResponseTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInput, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMax, double rAlpha) const;

    //! @brief ... calculate the (1-rAlpha) confidence interval on predictions of the linear regression
    //! @param rInput ... input point coordinates
    //! @param rOutputMean ... output of the regression model
    //! @param rOutputMin ... lower bound of the confidence interval
    //! @param rOutputMax ... upper bound of the confidence interval
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    void SolveConfidenceIntervalsPredictionTransformed(const FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rInput, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMean, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMin, NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rOutputMax, double rAlpha) const;

    //! @brief perform the regression significance test
    //! @param rTotalSumOfSquares ... total sum of squares
    //! @param rRegressionSumOfSquares ... regression sum of squares
    //! @param rErrorSumOfSquares ... error sum of squares
    //! @param rNumSupportPoints ... number of support points
    //! @param rNumRegressorVariables ... number of regressor variables
    //! @param rAlpha  ... the confidence level is defined as (1 - rAlpha)
    //! @return ... true if the regression is significant, false otherwise
    bool PerformRegressionSignificanceTest(double rTotalSumOfSquares, double rRegressionSumOfSquares, double rErrorSumOfSquares, int rNumSupportPoints, int rNumRegressorVariables, double rAlpha) const;

    //! @brief test the significance of a set of regressor variables/coefficients
    //! @param rTestCoefficients ... matrix of identifiers (number of identifiers, 1) defining the subset of regressor variables/coefficients (the intercept \f$\beta_0\f$ is not tested) for which the test is performed
    //! @param rAlpha ... the confidence level is defined as (1 - rAlpha)
    //! @param rTransformedFlag ... if true, the transformed support point coordinates are used, otherwise the original support point coordinates are used
    //! @return true if the contribution of the subset is significant, false otherwise
    bool PerformGeneralRegressionSignificanceTest(const NuTo::FullMatrix<int,Eigen::Dynamic,Eigen::Dynamic>& rTestCoefficients, double rAlpha, bool rTransformedFlag) const;

    //! @brief ... sum of squares of the error terms
	double Objective()const;

	//! @brief ... calculate the gradient of the objective function
	//! @param rGradient ... gradient vector (output)
	void Gradient(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& rGradient)const;

	//! @brief ... calculate the hessian of the objective function
	//! @param rHessian ... hessian matrix
    void Hessian(NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>&  rHessian)const;

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> mCoefficients; //!< ... regression coefficients \f$\beta_0, \beta_1, \ldots, \beta_p\f$
};

}
#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::MultipleLinearRegression)
#endif // SWIG
#endif // ENABLE_SERIALIZATION

