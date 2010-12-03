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

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>


#include "nuto/metamodel/MetamodelException.h"
#include "nuto/metamodel/MultipleLinearRegression.h"

// constructor
NuTo::MultipleLinearRegression::MultipleLinearRegression(): Metamodel()
{
}

// get regression coefficients
/*!
 * The regression variables are stored as vector
 *  \f[
 *  \boldsymbol{\beta} = \begin{bmatrix}
 *  \beta_0\\
 *  \beta_1\\
 *  \vdots\\
 *  \beta_p
 *  \end{bmatrix},
 *  \f]
 *  where \f$p\f$ is the number of regressor variables.
 */
void NuTo::MultipleLinearRegression::GetRegressionCoefficientsTransformed(NuTo::FullMatrix<double>& rRegressionCoefficients) const
{
	rRegressionCoefficients = this->mCoefficients;
}

// get covariance matrix of regression coefficients
/*!
 * The covariance matrix of the regression coefficients is defined as the inverse hessian matrix \f$\boldsymbol{H}\f$ multiplied by the mean square error \f$MS_e\f$
 * \f[
 *  COV(\boldsymbol{\beta}) = MS_e \boldsymbol{H}^{-1}.
 * \f]
 */
//! @sa NuTo::MultipleLinearRegression::Hessian(), NuTo::MultipleLinearRegression::GetMeanSquareErrorTransformed()
void NuTo::MultipleLinearRegression::GetRegressionCoefficientsCovarianceMatrixTransformed(NuTo::FullMatrix<double>& rCovarianceMatrix) const
{
	try
	{
		// calculate hessian matrix
		NuTo::FullMatrix<double> hessianMatrix;
		this->Hessian(hessianMatrix);

		// calculate inverse of hessian
		hessianMatrix.InverseCholeskyLapack(rCovarianceMatrix);

		// calculate covariance matrix
		rCovarianceMatrix *= this->GetMeanSquareErrorTransformed();
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetRegressionCoefficientsCovarianceMatrixTransformed] error calculating the covariance matrix of the regression coefficients.");
		throw myException;
	}
}

// get the (1-rAlpha) confidence interval on the regression coefficients
/*!
 * Assuming that the residuals (errors) \f$e_i\f$ of the regression model are normally and independetly distributed with zero mean and variance \f$\sigma_e^2\f$,
 * the confidence interval on the regression coefficients \f$\boldsymbol{\beta}\f$ are defined as
 * \f[
 *  \beta_i - t_{0.5\alpha,n-p} \sqrt{COV(\boldsymbol{\beta})_{ii}} \leq \beta_i \leq \beta_i + t_{0.5\alpha,n-p} \sqrt{COV(\boldsymbol{\beta})_{ii}}, \qquad i=0 \ldots p,
 * \f]
 * where \f$COV(\boldsymbol{\beta})\f$ is the covariance matrix of the regression coefficients, \f$n\f$ is the number of support points,
 * \f$p\f$ is the number regressor variables, and \f$t_{0.5\alpha,n-p}\f$ ist the \f$0.5\alpha\f$-quantile of the t-distribution with \f$n-p\f$ degrees of freedom.
 */
void NuTo::MultipleLinearRegression::GetRegressionCoefficientsConfidenceIntervalsTransformed(NuTo::FullMatrix<double>& rRegressionCoefficients, NuTo::FullMatrix<double>& rRegressionCoefficientsMin, NuTo::FullMatrix<double>& rRegressionCoefficientsMax, double rAlpha) const
{
	// get regression coefficients
	rRegressionCoefficients = this->mCoefficients;

	// get number of degrees of freedom
	int numDOF = this->mSupportPoints.GetNumSupportPoints() - this->mCoefficients.GetNumRows();
	if(numDOF < 1)
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionCoefficientsConfidenceIntervalsTransformed] number of support points must be larger than the number of regression coefficients.");
	}

	// calculate quantile of the t-distribution
	boost::math::students_t tDistribution(numDOF);
	double tValue = boost::math::quantile(boost::math::complement(tDistribution, 0.5 * rAlpha));

	try
	{
		// calculate covariance matrix
		NuTo::FullMatrix<double> covarianceMatrix;
		this->GetRegressionCoefficientsCovarianceMatrixTransformed(covarianceMatrix);

		// calculate confidence intervals
		rRegressionCoefficientsMin.Resize(this->mCoefficients.GetNumRows(),1);
		rRegressionCoefficientsMax.Resize(this->mCoefficients.GetNumRows(),1);
		for(int coefficient = 0; coefficient < this->mCoefficients.GetNumRows(); coefficient++)
		{
			double value = tValue * sqrt(covarianceMatrix.GetValue(coefficient, coefficient));
			rRegressionCoefficientsMin.SetValue(coefficient, 0, this->mCoefficients.GetValue(coefficient, 0) - value);
			rRegressionCoefficientsMax.SetValue(coefficient, 0, this->mCoefficients.GetValue(coefficient, 0) + value);
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetRegressionCoefficientsConfidenceIntervalsTransformed] error calculating confidence intervals on the regression coefficients.");
		throw myException;
	}
}

// calculate the (1-rAlpha) confidence interval on the mean response of the linear regression
void NuTo::MultipleLinearRegression::SolveConfidenceIntervalsMeanResponse(const FullMatrix<double>& rInput, NuTo::FullMatrix<double>& rOutputMean, NuTo::FullMatrix<double>& rOutputMin, NuTo::FullMatrix<double>& rOutputMax, double rAlpha) const
{
	try
	{
		//apply transformation of inputs
		NuTo::FullMatrix<double> rInputTransformed = rInput;
		mSupportPoints.TransformForwardInput(rInputTransformed);

		//solve the submodule
		SolveConfidenceIntervalsMeanResponseTransformed(rInputTransformed, rOutputMean, rOutputMin, rOutputMax, rAlpha);

		//apply transformation of outputs
		mSupportPoints.TransformForwardOutput(rOutputMean);
		mSupportPoints.TransformForwardOutput(rOutputMin);
		mSupportPoints.TransformForwardOutput(rOutputMax);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::SolveConfidenceIntervalsMeanResponse] error calculating confidence intervals on the mean response.");
		throw myException;
	}
}

// calculate the (1-rAlpha) confidence interval on the mean response of the linear regression (transformed data)
/*!
 * Assuming that the residuals (errors) \f$e_i\f$ of the regression model are normally and independetly distributed with zero mean and variance \f$\sigma_e^2\f$,
 * the confidence interval on the mean response \f$\hat{y}_i\f$ at point \f$\left[\hat{x}_{1i}, \hat{x}_{2i}, \ldots, \hat{x}_{pi}\right]\f$ is defined as
 * \f[
 *  \hat{y}_i - t_{0.5\alpha,n-p} \sqrt{\boldsymbol{\hat{x}}_i^T\boldsymbol{COV}(\boldsymbol{\beta}) \boldsymbol{\hat{x}}_i} \leq \hat{y}_i \leq  \hat{y}_i + t_{0.5\alpha,n-p} \sqrt{\boldsymbol{\hat{x}}_i^T\boldsymbol{COV}(\boldsymbol{\beta}) \boldsymbol{\hat{x}}_i}
 * \f]
 * with
 * \f{align*}{
 *  \hat{y}_i &= \beta_0 + \sum_{j=1}^p \beta_j \hat{x}_{ji}\\
 *  \boldsymbol{\hat{x}}_i &= \begin{bmatrix}
 *  1\\
 *  \hat{x}_{1i}\\
 *  \hat{x}_{2i}\\
 *  \vdots\\
 *  \hat{x}_{pi}
 *  \end{bmatrix},
 * \f}
 * where \f$COV(\boldsymbol{\beta})\f$ is the covariance matrix of the regression coefficients, \f$n\f$ is the number of support points,
 * \f$p\f$ is the number regressor variables, and \f$t_{0.5\alpha,n-p}\f$ ist the \f$0.5\alpha\f$-quantile of the t-distribution with \f$n-p\f$ degrees of freedom.
 */
void NuTo::MultipleLinearRegression::SolveConfidenceIntervalsMeanResponseTransformed(const FullMatrix<double>& rInput, NuTo::FullMatrix<double>& rOutputMean, NuTo::FullMatrix<double>& rOutputMin, NuTo::FullMatrix<double>& rOutputMax, double rAlpha) const
{
	// get number of degrees of freedom
	int numDOF = this->mSupportPoints.GetNumSupportPoints() - this->mCoefficients.GetNumRows();
	if(numDOF < 1)
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::GetSupportPointsConfidenceIntervalsTransformed] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate response
		this->SolveTransformed(rInput, rOutputMean);

		// calculate quantile of the t-distribution
		boost::math::students_t tDistribution(numDOF);
		double tValue = boost::math::quantile(boost::math::complement(tDistribution, 0.5 * rAlpha));

		// calculate covariance matrix
		NuTo::FullMatrix<double> covarianceMatrix;
		this->GetRegressionCoefficientsCovarianceMatrixTransformed(covarianceMatrix);

		// prepare output
		rOutputMin = rOutputMean;
		rOutputMax = rOutputMean;

		// calculate interval
		const double* curInputData = rInput.mEigenMatrix.data();
		for(int sample = 0; sample < rInput.GetNumColumns(); sample++)
		{
			const double* curCovarianceData = covarianceMatrix.mEigenMatrix.data();
			double sampleInterval = 0.0;
			for(int col = 0; col < covarianceMatrix.GetNumColumns(); col++)
			{
				double tmpValue = curCovarianceData[0];
				for(int row = 0; row < this->mSupportPoints.GetDimInput(); row++)
				{
					tmpValue += curCovarianceData[row + 1] * curInputData[row];
				}
				if(col == 0)
				{
					sampleInterval += tmpValue;
				}
				else
				{
					sampleInterval += tmpValue * curInputData[col - 1];
				}
				curCovarianceData += covarianceMatrix.GetNumRows();
			}
			sampleInterval = tValue * sqrt(sampleInterval);
			rOutputMin(0,sample) -= sampleInterval;
			rOutputMax(0,sample) += sampleInterval;

			curInputData += this->mSupportPoints.GetDimInput();
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::SolveConfidenceIntervalsMeanResponseTransformed] error calculating confidence intervals on the transformed mean of the response.");
		throw myException;
	}
}

// calculate the (1-rAlpha) confidence interval on predictions of the linear regression
void NuTo::MultipleLinearRegression::SolveConfidenceIntervalsPrediction(const FullMatrix<double>& rInput, NuTo::FullMatrix<double>& rOutputMean, NuTo::FullMatrix<double>& rOutputMin, NuTo::FullMatrix<double>& rOutputMax, double rAlpha) const
{
	try
	{
		//apply transformation of inputs
		NuTo::FullMatrix<double> rInputTransformed = rInput;
		mSupportPoints.TransformForwardInput(rInputTransformed);

		//solve the submodule
		SolveConfidenceIntervalsPredictionTransformed(rInputTransformed, rOutputMean, rOutputMin, rOutputMax, rAlpha);

		//apply transformation of outputs
		mSupportPoints.TransformForwardOutput(rOutputMean);
		mSupportPoints.TransformForwardOutput(rOutputMin);
		mSupportPoints.TransformForwardOutput(rOutputMax);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::SolveConfidenceIntervalsPrediction] error calculating confidence intervals on predictions of the response.");
		throw myException;
	}
}

// calculate the (1-rAlpha) confidence interval on predictions of the linear regression (transformed data)
/*!
 * Assuming that the residuals (errors) \f$e_i\f$ of the regression model are normally and independetly distributed with zero mean and variance \f$\sigma_e^2\f$,
 * the confidence interval on predictions \f$\hat{y}_i\f$ at point \f$\left[\hat{x}_{1i}, \hat{x}_{2i}, \ldots, \hat{x}_{pi}\right]\f$ is defined as
 * \f[
 *  \hat{y}_i - t_{0.5\alpha,n-p} \sqrt{\sigma_e^2 + \boldsymbol{\hat{x}}_i^T\boldsymbol{COV}(\boldsymbol{\beta}) \boldsymbol{\hat{x}}_i} \leq \hat{y}_i \leq  \hat{y}_i + t_{0.5\alpha,n-p} \sqrt{\sigma_e^2 + \boldsymbol{\hat{x}}_i^T\boldsymbol{COV}(\boldsymbol{\beta}) \boldsymbol{\hat{x}}_i}
 * \f]
 * with
 * \f{align*}{
 *  \hat{y}_i &= \beta_0 + \sum_{j=1}^p \beta_j \hat{x}_{ji}\\
 *  \boldsymbol{\hat{x}}_i &= \begin{bmatrix}
 *  1\\
 *  \hat{x}_{1i}\\
 *  \hat{x}_{2i}\\
 *  \vdots\\
 *  \hat{x}_{pi}
 *  \end{bmatrix},
 * \f}
 * where \f$COV(\boldsymbol{\beta})\f$ is the covariance matrix of the regression coefficients, \f$n\f$ is the number of support points,
 * \f$p\f$ is the number regressor variables, and \f$t_{0.5\alpha,n-p}\f$ ist the \f$0.5\alpha\f$-quantile of the t-distribution with \f$n-p\f$ degrees of freedom.
 * The mean-square-error of the support point residuals is used as estimator for the error variance \f$\sigma_e^2\f$.
 */
void NuTo::MultipleLinearRegression::SolveConfidenceIntervalsPredictionTransformed(const FullMatrix<double>& rInput, NuTo::FullMatrix<double>& rOutputMean, NuTo::FullMatrix<double>& rOutputMin, NuTo::FullMatrix<double>& rOutputMax, double rAlpha) const
{
	// get number of degrees of freedom
	int numDOF = this->mSupportPoints.GetNumSupportPoints() - this->mCoefficients.GetNumRows();
	if(numDOF < 1)
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::GetSupportPointsConfidenceIntervalsTransformed] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate response
		this->SolveTransformed(rInput, rOutputMean);

		// calculate quantile of the t-distribution
		boost::math::students_t tDistribution(numDOF);
		double tValue = boost::math::quantile(boost::math::complement(tDistribution, 0.5 * rAlpha));

		// calculate covariance matrix
		NuTo::FullMatrix<double> covarianceMatrix;
		this->GetRegressionCoefficientsCovarianceMatrixTransformed(covarianceMatrix);

		// calculate mean square error
		double mse = this->GetMeanSquareErrorTransformed();

		// prepare output
		rOutputMin = rOutputMean;
		rOutputMax = rOutputMean;

		// calculate interval
		const double* curInputData = rInput.mEigenMatrix.data();
		for(int sample = 0; sample < rInput.GetNumColumns(); sample++)
		{
			const double* curCovarianceData = covarianceMatrix.mEigenMatrix.data();
			double sampleInterval = 0.0;
			for(int col = 0; col < covarianceMatrix.GetNumColumns(); col++)
			{
				double tmpValue = curCovarianceData[0];
				for(int row = 0; row < this->mSupportPoints.GetDimInput(); row++)
				{
					tmpValue += curCovarianceData[row + 1] * curInputData[row];
				}
				if(col == 0)
				{
					sampleInterval += tmpValue;
				}
				else
				{
					sampleInterval += tmpValue * curInputData[col - 1];
				}
				curCovarianceData += covarianceMatrix.GetNumRows();
			}
			sampleInterval = tValue * sqrt(mse + sampleInterval);
			rOutputMin(0,sample) -= sampleInterval;
			rOutputMax(0,sample) += sampleInterval;

			curInputData += this->mSupportPoints.GetDimInput();
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::SolveConfidenceIntervalsPredictionTransformed] error calculating confidence intervals on predictions of the response (transformed).");
		throw myException;
	}
}

// calculate parameters
/*!
 * The regression coefficients are determined using the method of least squares for a given set of support points.
 * Due to the linear character of the regression model, the minimization problem can be directly solved for the unknown regression coefficients of.
 * The corresponding linear system of equations reads
 * \f[
 *  \boldsymbol{H} \boldsymbol{\beta} = - \left. \boldsymbol{G}\right|_{\beta_0=\beta_1=\ldots=\beta_p=0},
 * \f]
 * where \f$\boldsymbol{H}\f$ is the Hessian matrix of the objective function (constant due to the linear influence of the regression coefficients),
 * \f$\boldsymbol{G}\f$ is the gradient of the objective function (calculated for \f$\beta_0=\beta_1=\ldots=\beta_p=0\f$),
 * \f$p\f$ is the number of regressor variables, and \f$\boldsymbol{\beta}\f$ is the vector of regression coefficients.
 */
//! @sa NuTo::MultipleLinearRegression::Gradient(), NuTo::MultipleLinearRegression::Hessian()
void NuTo::MultipleLinearRegression::BuildDerived()
{
	// check dimension of output
	int dimOutput = this->mSupportPoints.GetDimOutput();
	if(dimOutput != 1)
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::BuildDerived] dimension of output must be 1.");
	}

	try
	{
		// resize parameter vector
		int dimInput = this->mSupportPoints.GetDimInput();
		int numParameters = dimInput + 1;
		this->mCoefficients.Resize(numParameters,1);
		this->mCoefficients.mEigenMatrix.setZero(numParameters,1);

		// get hessian matrix
		NuTo::FullMatrix<double> hessianMatrix;
		this->Hessian(hessianMatrix);

		// get gradient
		NuTo::FullMatrix<double> gradientVector;
		this->Gradient(gradientVector);
		gradientVector *= -1;

		// solve system of equations ( since hessian is constant -> only one step)
		hessianMatrix.SolveCholeskyLapack(gradientVector, this->mCoefficients);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::BuildDerived] error calculating regression coefficients.");
		throw myException;
	}
}

// calculate response for given input data (transformed system)
void NuTo::MultipleLinearRegression::SolveTransformed(const FullMatrix<double>& rInputCoordinates, NuTo::FullMatrix<double>& rOutputCoordinates)const
{
    int dimInput = this->mSupportPoints.GetDimInput();
	if(rInputCoordinates.GetNumRows() != dimInput)
	{
        throw MetamodelException("[NuTo::MultipleLinearRegression::SolveTransformed] Dimension of input (number of rows) is not identical with metamodel.");
	}
	int numCoefficients = this->mCoefficients.GetNumRows();
	if(numCoefficients != dimInput + 1)
	{
        throw MetamodelException("[NuTo::MultipleLinearRegression::SolveTransformed] invalid number of regression coefficients. Build model first.");
	}

	try
	{
		// prepare output
		int numSamples = rInputCoordinates.GetNumColumns();
		rOutputCoordinates.Resize(1, numSamples);

		// calculate output
		const double* coefficients = this->mCoefficients.mEigenMatrix.data();
		const double* inputData =  rInputCoordinates.mEigenMatrix.data();
		for(int sample = 0; sample < numSamples; sample++)
		{
			double curValue = coefficients[0];
			for(int parameter = 1; parameter < numCoefficients; parameter++)
			{
				curValue += coefficients[parameter] * inputData[0];
				inputData++;
			}
			rOutputCoordinates.SetValue(0, sample, curValue);
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::SolveTransformed] error calculating response.");
		throw myException;
	}
}

// calculate the residual between the support point data and the linear regression (original data)
/*!
 * The residual \f$e\f$ is defined as the difference between the support point outputs and the corresponding outputs of the linear regression model
 * \f[
 *   e_i = y_i - \beta_0 - \sum_{j=1}^p \beta_j x_{ji},\qquad i = 1 \ldots n
 * \f]
 * where \f$n\f$ is the number of support points, \f$p\f$ is the number of regression variables, \f$\beta_0, \beta_1, \ldots, \beta_p\f$
 *  are the regression coefficients, and
 *  \f[\begin{bmatrix}
 *  y_1 & x_{11} & x_{21} & \ldots & x_{p1}\\
 *  y_2 & x_{12} & x_{22} & \ldots & x_{p2}\\
 *  \vdots & \vdots & \vdots & \ddots & \vdots\\
 *  y_n & x_{1n} & x_{2n} & \ldots & x_{pn}
 *  \end{bmatrix}\f]
 *  are the support points.
 */
void NuTo::MultipleLinearRegression::GetSupportPointsResidual(NuTo::FullMatrix<double>& rResidualVector) const
{
	// check output dimension
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetSupportPointsResidual] number of output variables must be one.");
	}

	try
	{
		// calculate output of the regression model
		NuTo::FullMatrix<double> regressionOutputVector;
		this->SolveTransformed(this->mSupportPoints.GetTransformedSupportPointsInput(), regressionOutputVector);

		// apply transformation of outputs
		this->mSupportPoints.TransformForwardOutput(regressionOutputVector);

		// substract original output
		rResidualVector = this->mSupportPoints.GetOrigSupportPointsOutput() - regressionOutputVector;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetSupportPointsResidual] error calculating residuals.");
		throw myException;
	}
}

// calculate the residual between the support point data and the linear regression (transformed data)
//! @sa NuTo::MultipleLinearRegression::GetSupportPointsResidual
void NuTo::MultipleLinearRegression::GetSupportPointsResidualTransformed(NuTo::FullMatrix<double>& rResidualVector) const
{
	// check output dimension
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetSupportPointsResidualTransformed] number of output variables must be one.");
	}

	try
	{
		// calculate output of the regression model
		NuTo::FullMatrix<double> regressionOutputVector;
		this->SolveTransformed(this->mSupportPoints.GetTransformedSupportPointsInput(), regressionOutputVector);

		// substract original output
		rResidualVector = this->mSupportPoints.GetTransformedSupportPointsOutput() - regressionOutputVector;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetSupportPointsResidualTransformed] error calculating transformed residuals.");
		throw myException;
	}
}

// calculate the mean square error using transformed support points
/*!
 * The mean square error \f$MS_e\f$, which is an estimator of the error variance \f$\sigma_e^2\f$, is defined as
 * \f[
 * \sigma_e^2 = MS_e = \dfrac{1}{n-p-1}\sum_{i=1}^n e_i^2,
 * \f]
 * where \f$e_i\f$ are the residuals (difference between the support point outputs and the corresponding regression outputs),
 * \f$n\f$ is the number of support points, and \f$p\f$ is the number of regressor variables.
 */
double NuTo::MultipleLinearRegression::GetMeanSquareErrorTransformed() const
{
	// check data
	if(this->mSupportPoints.GetNumSupportPoints() <= this->mCoefficients.GetNumRows())
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::GetMeanSquareErrorTransformed] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate mean square error
		return this->GetErrorSumOfSquaresTransformed()/static_cast<double>(this->mSupportPoints.GetNumSupportPoints() - this->mCoefficients.GetNumRows());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetMeanSquareErrorTransformed] error calculating mean square error.");
		throw myException;
	}
}

// calculate the mean square error using original support points
/*!
 * The mean square error \f$MS_e\f$, which is an estimator of the error variance \f$\sigma_e^2\f$, is defined as
 * \f[
 * \sigma_e^2 = MS_e = \dfrac{1}{n-p-1}\sum_{i=1}^n e_i^2,
 * \f]
 * where \f$e_i\f$ are the residuals (difference between the support point outputs and the corresponding regression outputs),
 * \f$n\f$ is the number of support points, and \f$p\f$ is the number of regressor variables.
 */
//! @sa NuTo::MultipleLinearRegression::GetSupportPointsResidual
double NuTo::MultipleLinearRegression::GetMeanSquareError() const
{
	// check data
	if(this->mSupportPoints.GetNumSupportPoints() <= this->mCoefficients.GetNumRows())
	{
		throw MetamodelException("[NuTo::MultipleLinearRegression::GetMeanSquareError] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate mean square error
		return this->GetErrorSumOfSquares()/static_cast<double>(this->mSupportPoints.GetNumSupportPoints() - this->mCoefficients.GetNumRows());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetMeanSquareError] error calculating mean square error.");
		throw myException;
	}
}

// calculate coefficient of determination using original support point data
/*!
 * The coefficient of determination is defined \f$R^2\f$ as
 * \f[
 *  R^2 = \dfrac{SS_r}{S_{yy}} = 1 - \dfrac{SS_e}{S_{yy}},
 * \f]
 * where \f$S_{yy}\f$ is the total sum of squares, \f$SS_e\f$ is the error sum of squares, and \f$SS_r\f$ is the regression sum of squares.
 */
double NuTo::MultipleLinearRegression::GetCoefficientOfDetermination() const
{
	try
	{
		return this->GetRegressionSumOfSquares()/this->GetTotalSumOfSquares();
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetCoefficientOfDetermination] error calculating coefficient of determination.");
		throw myException;
	}
}

// calculate coefficient of determination using transformed support point data
/*!
 * The coefficient of determination is defined \f$R^2\f$ as
 * \f[
 *  R^2 = \dfrac{SS_r}{S_{yy}} = 1 - \dfrac{SS_e}{S_{yy}},
 * \f]
 * where \f$S_{yy}\f$ is the total sum of squares, \f$SS_e\f$ is the error sum of squares, and \f$SS_r\f$ is the regression sum of squares.
 */
double NuTo::MultipleLinearRegression::GetCoefficientOfDeterminationTransformed() const
{
	try
	{
		return this->GetRegressionSumOfSquaresTransformed()/this->GetTotalSumOfSquaresTransformed();
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetCoefficientOfDeterminationTransformed] error calculating coefficient of determination.");
		throw myException;
	}
}

// calculate adjusted coefficient of determination
/*!
 * The adjusted coefficient of determination \f$R^2_{adj}\f$ is defined as
 * \f[
 *  R^2_{adj} = 1.0 - (1.0 - R^2)\dfrac{n - 1}{n - p - 1},
 * \f]
 * where \f$R^2\f$ is the coefficient of determination, \f$n\f$ is the number of support points,
 * and \f$p\f$ is the number of regressor variables.
 */
double NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDetermination() const
{
	// check variables
	double numSamples = this->mSupportPoints.GetNumSupportPoints();
	double numUnknowns = this->mCoefficients.GetNumRows();
	if(numSamples - numUnknowns < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDetermination] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate adjusted coefficient of determination
		return 1.0 - (numSamples - 1.0)/(numSamples - numUnknowns)*(1.0 - this->GetCoefficientOfDetermination());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDetermination] error adjusted calculating coefficient of determination.");
		throw myException;
	}
}

// calculate adjusted coefficient of determination
/*!
 * The adjusted coefficient of determination \f$R^2_{adj}\f$ is defined as
 * \f[
 *  R^2_{adj} = 1.0 - (1.0 - R^2)\dfrac{n - 1}{n - p - 1},
 * \f]
 * where \f$R^2\f$ is the coefficient of determination, \f$n\f$ is the number of support points,
 * and \f$p\f$ is the number of regressor variables.
 */
double NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDeterminationTransformed() const
{
	// check variables
	double numSamples = this->mSupportPoints.GetNumSupportPoints();
	double numUnknowns = this->mCoefficients.GetNumRows();
	if(numSamples - numUnknowns < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDeterminationTransformed] number of support points must be larger than the number of regression coefficients.");
	}

	try
	{
		// calculate adjusted coefficient of determination
		return 1.0 - (numSamples - 1.0)/(numSamples - numUnknowns)*(1.0 - this->GetCoefficientOfDeterminationTransformed());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetAdjustedCoefficientOfDeterminationTransformed] error adjusted calculating coefficient of determination.");
		throw myException;
	}
}

// return type ID
std::string NuTo::MultipleLinearRegression::GetTypeId()const
{
    return std::string("MultipleLinearRegression");
}

// calculate objective
/*!
 *  Using the method of least squares, the objective \f$L\f$ is defined as the sum of squares of the error terms \f$e_j\f$
 *  \f[L=\dfrac{1}{2}\sum_{j=1}^n \left[e_j\right]^2=\dfrac{1}{2}\sum_{j=1}^n\left[y_j - \beta_0 - \sum_{i=1}^p \beta_i x_{ij}\right]^2,\f]
 *  where \f$n\f$ is the number of support points, \f$p\f$ is the number of regression variables, \f$\beta_0, \beta_1, \ldots, \beta_p\f$
 *  are the regression coefficients, and
 *  \f[\begin{bmatrix}
 *  y_1 & x_{11} & x_{21} & \ldots & x_{p1}\\
 *  y_2 & x_{12} & x_{22} & \ldots & x_{p2}\\
 *  \vdots & \vdots & \vdots & \ddots & \vdots\\
 *  y_n & x_{1n} & x_{2n} & \ldots & x_{pn}
 *  \end{bmatrix}\f]
 *  are the support points.
 */
double NuTo::MultipleLinearRegression::Objective()const
{
	try
	{
		// get support point output
		assert(this->mSupportPoints.GetDimOutput() == 1);
		const double* outputData = this->mSupportPoints.GetTransformedSupportPointsOutput().mEigenMatrix.data();

		// get support point input
		int dimInput = this->mSupportPoints.GetDimInput();
		int numCoefficients = dimInput + 1;
		assert(this->mCoefficients.size() == numCoefficients);
		const double* inputData = this->mSupportPoints.GetTransformedSupportPointsInput().mEigenMatrix.data();

		// get number of samples
		int numSamples = this->mSupportPoints.GetNumSupportPoints();

		// calculate objective
		double objective = 0.0;
		const double* curInputData = inputData;
		const double* coefficients = this->mCoefficients.mEigenMatrix.data();
		for(int sample = 0; sample < numSamples; sample++)
		{
			double term = outputData[sample] - coefficients[0];
			for(int variable = 1; variable < numCoefficients; variable++)
			{
				term -= coefficients[variable] * curInputData[variable - 1];
			}
			objective += term * term;
			curInputData += dimInput;
		}
		return 0.5 * objective;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::Objective] error calculating objective.");
		throw myException;
	}
}

// calculate gradient of the objective function
/*!
 *  The gradient of the objective function \f$L\f$ is defined as the partial derivatives of the sum of the error terms
 *  with respect to the regression coefficients
 *  \f{align*}{
 *  \dfrac{\partial L}{\partial \beta_0} &= \sum_{j=1}^n \left[\beta_0 - y_j + \sum_{i=1}^p \beta_i x_{ij} \right]&
 *  &=n \beta_0 - \sum_{j=1}^n y_j + \sum_{j=1}^n \sum_{i=1}^p \beta_i x_{ij} & \\
 *  \dfrac{\partial L}{\partial \beta_k} &= \sum_{j=1}^n \left[\beta_0 - y_j + \sum_{i=1}^p \beta_i x_{ij} \right] x_{kj}&
 *  &= \beta_0 \sum_{j=1}^n x_{kj} - \sum_{j=1}^n y_j x_{kj} + \sum_{i=1}^p \beta_i \sum_{j=1}^n x_{ij} x_{kj} &\qquad k=1 \ldots p
 *  \f}
 *  where \f$n\f$ is the number of support points, \f$p\f$ is the number of regression variables, \f$\beta_0, \beta_1, \ldots, \beta_p\f$
 *  are the regression coefficients, and
 *  \f[\begin{bmatrix}
 *  y_1 & x_{11} & x_{21} & \ldots & x_{p1}\\
 *  y_1 & x_{12} & x_{22} & \ldots & x_{p2}\\
 *  \vdots & \vdots & \vdots & \ddots & \vdots\\
 *  y_n & x_{1n} & x_{2n} & \ldots & x_{pn}
 *  \end{bmatrix}\f]
 *  are the support points. The gradient is stored as full matrix with number of regression coefficients rows and one column:
 *  \f[
 *  \boldsymbol{G} = \begin{bmatrix}
 *  \dfrac{\partial L}{\partial \beta_0} \\
 *  \dfrac{\partial L}{\partial \beta_1} \\
 *  \cdots \\
 *  \dfrac{\partial L}{\partial \beta_k}
 *  \end{bmatrix}.
 *  \f]
 */
void NuTo::MultipleLinearRegression::Gradient(NuTo::FullMatrix<double>& rGradient)const
{
	try
	{
		// get support point output
		assert(this->mSupportPoints.GetDimOutput() == 1);
		const double* outputData = this->mSupportPoints.GetTransformedSupportPointsOutput().mEigenMatrix.data();

		// get support point input
		int dimInput = this->mSupportPoints.GetDimInput();
		int numCoefficients = dimInput + 1;
		assert(this->mCoefficients.size() == numCoefficients);
		const double* inputData = this->mSupportPoints.GetTransformedSupportPointsInput().mEigenMatrix.data();

		// get number of samples
		int numSamples = this->mSupportPoints.GetNumSupportPoints();

		// prepare gradient
		rGradient.Resize(numCoefficients, 1);

		// get coefficients
		const double* coefficients = this->mCoefficients.mEigenMatrix.data();

		// calculate derivative with respect to beta_0
		double curDerivative = numSamples * coefficients[0];
		for(int sample = 0; sample < numSamples; sample++)
		{
			curDerivative -= outputData[sample];
		}
		for(int variable = 1; variable < numCoefficients; variable++)
		{
			if(coefficients[variable] != 0.0)
			{
				const double *curInputData = inputData + variable - 1;
				for(int sample = 0; sample < numSamples; sample++)
				{
					curDerivative += coefficients[variable] *  curInputData[0];
					curInputData += dimInput;
				}
			}
		}
		rGradient.SetValue(0,0, curDerivative);

		// calculate derivatives with respect to beta_k (k=1..p)
		for(int row = 1; row < numCoefficients; row++)
		{
			double curDerivative = 0.0;
			if(coefficients[0] == 0.0)
			{
				const double* curInputData = inputData + row - 1;
				for(int sample = 0; sample < numSamples; sample++)
				{
					curDerivative -= outputData[sample] * curInputData[0];
					curInputData += dimInput;
				}
			}
			else
			{
				const double* curInputData = inputData + row - 1;
				for(int sample = 0; sample < numSamples; sample++)
				{
					curDerivative += (coefficients[0] - outputData[sample]) * curInputData[0];
					curInputData += dimInput;
				}
			}
			for(int variable = 1; variable < numCoefficients; variable++)
			{
				if(coefficients[variable] != 0.0)
				{
					const double* curInputDataRow = inputData + row - 1;
					const double* curInputDataVariable = inputData + variable - 1;
					for(int sample = 0; sample < numSamples; sample++)
					{
						curDerivative += coefficients[variable] * curInputDataRow[0] * curInputDataVariable[0];
						curInputDataRow += dimInput;
						curInputDataVariable += dimInput;
					}
				}
			}
			rGradient.SetValue(row,0, curDerivative);
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::Gradient] error calculating gradient.");
		throw myException;
	}
}

// calculate hessian matrix
/*!
 * The Hessian matrix of the objective function \f$L\f$ is defined as the second drivative of the sum of squares of error terms
 *  with respect to the regression coefficients
 *  \f{align*}{
 *   \dfrac{\partial^2 L}{\partial \beta_0^2} &= \dfrac{1}{n} &\\
 *   \dfrac{\partial^2 L}{\partial \beta_0 \partial \beta_l} &= \sum_{j=1}^n x_{lj} & \qquad l = 1, 2, \ldots, p\\
 *   \dfrac{\partial^2 L}{\partial \beta_l \partial \beta_k} &= \sum_{j=1}^n x_{lj} x_{kj} & \qquad k = 1, 2, \ldots, p,
 *  \f}
 *  where \f$n\f$ is the number of support points, \f$p\f$ is the number of regression variables, \f$\beta_0, \beta_1, \ldots, \beta_p\f$
 *  are the regression coefficients, and
 *  \f[\begin{bmatrix}
 *  y_1 & x_{11} & x_{21} & \ldots & x_{p1}\\
 *  y_2 & x_{12} & x_{22} & \ldots & x_{p2}\\
 *  \vdots & \vdots & \vdots & \ddots & \vdots\\
 *  y_n & x_{1n} & x_{2n} & \ldots & x_{pn}
 *  \end{bmatrix}\f]
 *  are the support points. The symmetric hessian matrix is stored as full matrix (dimension: number of regression coefficients):
 *  \f[\boldsymbol{H} =
 *  \begin{bmatrix}
 *  \dfrac{\partial^2 L}{\partial \beta_0^2} & \dfrac{\partial^2 L}{\partial \beta_0 \partial \beta_1} & \cdots & \dfrac{\partial^2 L}{\partial \beta_0 \partial \beta_p} \\
 *  \dfrac{\partial^2 L}{\partial \beta_0 \partial \beta_1} & \dfrac{\partial^2 L}{\partial \beta_1^2} & \cdots & \dfrac{\partial^2 L}{\partial \beta_1 \partial \beta_p} \\
 *  \vdots & \vdots & \ddots & \vdots\\
 *  \dfrac{\partial^2 L}{\partial \beta_0 \partial \beta_p} & \dfrac{\partial^2 L}{\partial \beta_1 \partial \beta_p} & \cdots & \dfrac{\partial^2 L}{\partial \beta_p^2}
 *  \end{bmatrix}.
 *  \f]
 */
void NuTo::MultipleLinearRegression::Hessian(NuTo::FullMatrix<double>&  rHessian)const
{
	try
	{
		// check output
		assert(this->mSupportPoints.GetDimOutput() == 1);

		// get support point input
		int dimInput = this->mSupportPoints.GetDimInput();
		int numCoefficients = dimInput + 1;
		assert(this->mCoefficients.size() == numCoefficients);
		const double* inputData = this->mSupportPoints.GetTransformedSupportPointsInput().mEigenMatrix.data();

		// get number of samples
		int numSamples = this->mSupportPoints.GetNumSupportPoints();

		// prepare hessian
		rHessian.Resize(numCoefficients,numCoefficients);

		// calculate Hessian
		rHessian.SetValue(0, 0, numSamples);

		for(int coefficient = 1; coefficient < numCoefficients; coefficient++)
		{
			double curDerivative = 0.0;
			const double* curInputData = inputData + coefficient - 1;
			for(int sample = 0; sample < numSamples; sample++)
			{
				curDerivative += curInputData[0];
				curInputData += dimInput;
			}
			rHessian.SetValue(coefficient, 0, curDerivative);
			rHessian.SetValue(0, coefficient, curDerivative);
		}

		for(int row = 1; row < numCoefficients; row++)
		{
			for(int col = row; col < numCoefficients; col++)
			{
				double curDerivative = 0.0;
				const double* curRowInputData = inputData + row - 1;
				const double* curColInputData = inputData + col - 1;
				for(int sample = 0; sample < numSamples; sample++)
				{
					curDerivative += curRowInputData[0] * curColInputData[0];
					curRowInputData += dimInput;
					curColInputData += dimInput;
				}
				rHessian.SetValue(row, col, curDerivative);
				if(row != col)
				{
					rHessian.SetValue(col, row, curDerivative);
				}
			}
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::Hessian] error calculating hessian.");
		throw myException;
	}
}

// calculate total sum of squares
/*!
 * The total sum of squares \f$S_{yy}\f$ is defined as
 * \f[
 *  S_{yy} = \sum_{i=1}^n \left[y_i - \bar{y}\right]^2 = \sum_{i=1}^n y_i^2 - \dfrac{1}{n} \left(\sum_{i=1}^n y_i\right)^2,
 * \f]
 * where \f$y_i\f$ is the output of support point \f$i\f$, \f$\bar{y}\f$ is the corresponding mean value, and \f$n\f$ is the number of support points.
 */
double NuTo::MultipleLinearRegression::GetTotalSumOfSquares() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquares] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquares] invalid number of regression coefficients (build model first).");
	}
	if(this->mSupportPoints.GetNumSupportPoints() < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquares] invalid number of support points.");
	}

	try
	{
		// calculate total sum of squares
		double sum = 0.0;
		double sumSquare = 0.0;
		const double* outputData = this->mSupportPoints.GetOrigSupportPointsOutput().mEigenMatrix.data();
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sum += outputData[sample];
			sumSquare += outputData[sample] * outputData[sample];
		}
		return sumSquare - 1.0/static_cast<double>(this->mSupportPoints.GetNumSupportPoints()) * sum * sum;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetTotalSumOfSquares] error calculating total sum of squares.");
		throw myException;
	}
}

// calculate total sum of squares (transformed support points)
/*!
 * The total sum of squares \f$S_{yy}\f$ is defined as
 * \f[
 *  S_{yy} = \sum_{i=1}^n \left[y_i - \bar{y}\right]^2 = \sum_{i=1}^n y_i^2 - \dfrac{1}{n} \left(\sum_{i=1}^n y_i\right)^2,
 * \f]
 * where \f$y_i\f$ is the output of support point \f$i\f$, \f$\bar{y}\f$ is the corresponding mean value, and \f$n\f$ is the number of support points.
 */
double NuTo::MultipleLinearRegression::GetTotalSumOfSquaresTransformed() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquaresTransformed] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquaresTransformed] invalid number of regression coefficients (build model first).");
	}
	if(this->mSupportPoints.GetNumSupportPoints() < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetTotalSumOfSquaresTransformed] invalid number of support points.");
	}

	// calculate total sum of squares
	try
	{
		double sum = 0.0;
		double sumSquare = 0.0;
		const double* outputData = this->mSupportPoints.GetTransformedSupportPointsOutput().mEigenMatrix.data();
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sum += outputData[sample];
			sumSquare += outputData[sample] * outputData[sample];
		}
		return sumSquare - 1.0/static_cast<double>(this->mSupportPoints.GetNumSupportPoints()) * sum * sum;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetTotalSumOfSquaresTransformed] error calculating total sum of squares.");
		throw myException;
	}
}

// calculate error (residual) sum of squares
/*!
 * The error (residual) sum of squares \f$SS_e\f$ is defined as
 * \f[
 * SS_e = \sum_{i=1}^n e_i^2 = \sum_{i=1}^n \left[y_i - \beta_0 - \sum_{j=1}^p \beta_j x_{ji} \right]^2 = \sum_{i=1}^n y_i^2 - \sum_{i=1}^n \hat{y}_i y_i,
 * \f]
 * where \f$e_i\f$ is the rsidual, \f$n\f$ is the number of support points, \f$p\f$ is the number of regressor variables,
 * \f$\left[y_i, x_{0i}, x_{1i}, \ldots, x_{pi}\right]\f$ are the coordinates of the support point \f$i\f$.
 */
double NuTo::MultipleLinearRegression::GetErrorSumOfSquares() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetErrorSumOfSquares] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetErrorSumOfSquares] invalid number of regression coefficients (build model first).");
	}

	try
	{
		// calculate residuals
		NuTo::FullMatrix<double> residualVector;
		this->GetSupportPointsResidual(residualVector);

		// calculate sum of squares of residuals
		double sse = 0.0;
		const double* data = residualVector.mEigenMatrix.data();
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sse += data[sample] * data[sample];
		}
		return sse;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetErrorSumOfSquares] error calculating error sum of squares.");
		throw myException;
	}
}

// calculate error (residual) sum of squares (transformed support points)
/*!
 * The error (residual) sum of squares \f$SS_e\f$ is defined as
 * \f[
 * SS_e = \sum_{i=1}^n e_i^2 = \sum_{i=1}^n \left[y_i - \beta_0 - \sum_{j=1}^p \beta_j x_{ji} \right]^2 = \sum_{i=1}^n y_i^2 - \sum_{i=1}^n \hat{y}_i y_i,
 * \f]
 * where \f$e_i\f$ is the rsidual, \f$n\f$ is the number of support points, \f$p\f$ is the number of regressor variables,
 * \f$\left[y_i, x_{0i}, x_{1i}, \ldots, x_{pi}\right]\f$ are the coordinates of the support point \f$i\f$.
 */
double NuTo::MultipleLinearRegression::GetErrorSumOfSquaresTransformed() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetErrorSumOfSquaresTransformed] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetErrorSumOfSquaresTransformed] invalid number of regression coefficients (build model first).");
	}

	try
	{
		// calculate residuals
		NuTo::FullMatrix<double> residualVector;
		this->GetSupportPointsResidualTransformed(residualVector);

		// calculate sum of squares of residuals
		double sse = 0.0;
		const double* data = residualVector.mEigenMatrix.data();
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sse += data[sample] * data[sample];
		}
		return sse;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetErrorSumOfSquaresTransformed] error calculating error sum of squares.");
		throw myException;
	}
}

// get regression sum of squares (transformed support points)
/*!
 * The regression sum of squares \f$SS_r\f$ is obtained from the split of the total sum of squares \f$S_{yy}\f$ into a sum of squares due to regression
 * \f$SS_r\f$ and the sum of squares due to error
 * \f[
 * SS_r = S_{yy} - SS_e = \sum_{i=1}^n \hat{y}_i y_i - \dfrac{1}{n}\left[\sum_{i=1}^n y_i \right]^2.
 * \f]
 */
double NuTo::MultipleLinearRegression::GetRegressionSumOfSquaresTransformed() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquaresTransformed] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquaresTransformed] invalid number of regression coefficients (build model first).");
	}
	if(this->mSupportPoints.GetNumSupportPoints() < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquaresTransformed] invalid number of support points.");
	}

	try
	{
		// calculate regression for support points
		NuTo::FullMatrix<double> outputRegression;
		this->SolveTransformed(this->mSupportPoints.GetTransformedSupportPointsInput(), outputRegression);
		const double* outputRegressionData = outputRegression.mEigenMatrix.data();

		// get original support point data
		const double* outputSupportPointsData = this->mSupportPoints.GetTransformedSupportPointsOutput().mEigenMatrix.data();

		// calculate sums
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sum1 += outputSupportPointsData[sample] * outputRegressionData[sample];
			sum2 += outputSupportPointsData[sample];
		}

		// calculate regression sum of squares
		return sum1 - sum2 * sum2 / static_cast<double>(this->mSupportPoints.GetNumSupportPoints());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquaresTransformed] error calculating regression sum of squares.");
		throw myException;
	}

}

// get regression sum of squares
/*!
 * The regression sum of squares \f$SS_r\f$ is obtained from the split of the total sum of squares \f$S_{yy}\f$ into a sum of squares due to regression
 * \f$SS_r\f$ and the sum of squares due to error
 * \f[
 * SS_r = S_{yy} - SS_e = \sum_{i=1}^n \hat{y}_i y_i - \dfrac{1}{n}\left[\sum_{i=1}^n y_i \right]^2.
 * \f]
 */
double NuTo::MultipleLinearRegression::GetRegressionSumOfSquares() const
{
	// check input
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquares] dimension of output must be 1.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquares] invalid number of regression coefficients (build model first).");
	}
	if(this->mSupportPoints.GetNumSupportPoints() < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquares] invalid number of support points.");
	}

	try
	{
		// calculate regression for support points
		NuTo::FullMatrix<double> outputRegression;
		this->Solve(this->mSupportPoints.GetOrigSupportPointsInput(), outputRegression);
		const double* outputRegressionData = outputRegression.mEigenMatrix.data();

		// get original support point data
		const double* outputSupportPointsData = this->mSupportPoints.GetOrigSupportPointsOutput().mEigenMatrix.data();

		// calculate sums
		double sum1 = 0.0;
		double sum2 = 0.0;
		for(int sample = 0; sample < this->mSupportPoints.GetNumSupportPoints(); sample++)
		{
			sum1 += outputSupportPointsData[sample] * outputRegressionData[sample];
			sum2 += outputSupportPointsData[sample];
		}

		// calculate regression sum of squares
		return sum1 - sum2 * sum2 / static_cast<double>(this->mSupportPoints.GetNumSupportPoints());
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::GetRegressionSumOfSquares] error calculating regression sum of squares.");
		throw myException;
	}
}

// test for the significance of the regression
/*!
 * This hypothesis test determines whether there is a linear relationship between the dependent variable \f$y\f$ and the regressor variables \f$x_1, x_2, \ldots, x_p\f$.
 * The corresponding hypotheses are:
 * \f{align*}{
 * H_0: &\quad \beta_1 = \beta_2 = \ldots = \beta_p = 0.0,\\
 * H_1: &\quad \beta_i \neq 0.0 \text{ for at least one }i.
 * \f}
 * The null hypothesis \f$H_0\f$ is rejected (at least one regressor variable has a linear influence on the dependent variable) if
 * \f[
 *  F_0 = \dfrac{MS_r}{MS_e} = \dfrac{SS_r (n-p-1)}{SS_e p} > F_{\alpha,p,n-p-1},
 * \f]
 * where \f$SS_r\f$ is the regression sum of squares, \f$SS_e\f$ is the error sum of squares, \f$n\f$ is the number of support points,
 * \f$p\f$ is the number of regressor variables, \f$F_{\alpha,p,n-p-1}\f$ is the F-distribution.
 */
bool NuTo::MultipleLinearRegression::TestRegressionSignificance(double rAlpha) const
{
	try
	{
		// calculate sum of squares
		double Syy = this->GetTotalSumOfSquares();
		double SSe = this->GetErrorSumOfSquares();
		double SSr = this->GetRegressionSumOfSquares();

		// perform significance test
		return this->PerformRegressionSignificanceTest(Syy, SSr, SSe, this->mSupportPoints.GetNumSupportPoints(), this->mSupportPoints.GetDimInput(), rAlpha);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::TestRegressionSignificance] error testing significance of regression.");
		throw myException;
	}
}

// test for the significance of the regression (transformed support points)
/*!
 * This hypothesis test determines whether there is a linear relationship between the dependent variable \f$y\f$ and the regressor variables \f$x_1, x_2, \ldots, x_p\f$.
 * The corresponding hypotheses are:
 * \f{align*}{
 * H_0: &\quad \beta_1 = \beta_2 = \ldots = \beta_p = 0.0,\\
 * H_1: &\quad \beta_i \neq 0.0 \text{ for at least one }i.
 * \f}
 * The null hypothesis \f$H_0\f$ is rejected (at least one regressor variable has a linear influence on the dependent variable) if
 * \f[
 *  F_0 = \dfrac{MS_r}{MS_e} = \dfrac{SS_r (n-p-1)}{SS_e p} > F_{\alpha,p,n-p-1},
 * \f]
 * where \f$SS_r\f$ is the regression sum of squares, \f$SS_e\f$ is the error sum of squares, \f$n\f$ is the number of support points,
 * \f$p\f$ is the number of regressor variables, \f$F_{\alpha,p,n-p-1}\f$ is the F-distribution.
 */
bool NuTo::MultipleLinearRegression::TestRegressionSignificanceTransformed(double rAlpha) const
{
	try
	{
		// calculate sum of squares
		double Syy = this->GetTotalSumOfSquaresTransformed();
		double SSe = this->GetErrorSumOfSquaresTransformed();
		double SSr = this->GetRegressionSumOfSquaresTransformed();

		// perform significance test
		return this->PerformRegressionSignificanceTest(Syy, SSr, SSe, this->mSupportPoints.GetNumSupportPoints(), this->mSupportPoints.GetDimInput(), rAlpha);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::TestRegressionSignificanceTransformed] error testing significance of regression.");
		throw myException;
	}
}

// test for the significance of the individual regression coefficients
/*!
 * The hypotheses for testing the significance of the regression coefficient \f$\beta_i\f$ are
 * \f{align*}{
 * H_0: &\quad \beta_i = 0.0\\
 * H_1: &\quad \beta_i \neq 0.0.
 * \f}
 * The test statistic for this hypothesis is
 * \f[
 *  t_0 = \dfrac{\beta_i}{\sqrt{COV(\boldsymbol{\beta})_{ii}}},
 * \f]
 * where \f$\boldsymbol{COV(\boldsymbol{\beta})}\f$ is the covariance matrix of the regression coefficients \f$\boldsymbol{\beta}\f$.
 * The null hypothesis \f$H_0\f$ is rejected if
 * \f[
 *  \left| t_0 \right| > t_{0.5\alpha,n-p-1},
 * \f]
 * where \f$n\f$ is the number of support points, \f$p\f$ is the number of regressor variables, and \f$t\f$ is the distribution.
 */
void NuTo::MultipleLinearRegression::TestRegressionCoefficientsSignificanceTransformed(std::vector<bool>& rTestResult, double rAlpha) const
{
	try
	{
		// calculate covariance matrix
		NuTo::FullMatrix<double> covarianceMatrix;
		this->GetRegressionCoefficientsCovarianceMatrixTransformed(covarianceMatrix);

		// calculate rAlpha/2 quantile of t-distribution
		boost::math::students_t tDistribution(this->mSupportPoints.GetNumSupportPoints()-this->mSupportPoints.GetDimInput()-1);
		double t = boost::math::quantile(boost::math::complement(tDistribution, 0.5 * rAlpha));

		// prepare output
		if(this->mVerboseLevel > 0)
		{
			std::cout << std::endl;
			std::cout << "Multiple linear regression: regression coefficient significance test" << std::endl;
			std::cout << "  null hypothesis H0: beta_i = 0.0" << std::endl;
			std::cout << "  alternative hypothesis H1: beta_i != 0.0" << std::endl;
		}
		// perform tests
		rTestResult.resize(this->mCoefficients.GetNumRows());
		for(int coefficient = 0; coefficient < this->mCoefficients.GetNumRows(); coefficient++)
		{
			double abs_t0 = fabs(this->mCoefficients(coefficient,0)/sqrt(covarianceMatrix(coefficient,coefficient)));
			if(abs_t0 > t)
			{
				rTestResult[coefficient] = true;
			}
			else
			{
				rTestResult[coefficient] = false;
			}
			if(this->mVerboseLevel > 0)
			{
				std::cout << "  statistic beta_" << coefficient << ": |t0| = " << abs_t0;
				if(abs_t0 > t)
				{
					std::cout << " > t = " << t << " --> The null hypothesis H0 is rejected." << std::endl;
				}
				else
				{
					if(abs_t0 < t)
					{
						std::cout << " < ";
					}
					else
					{
						std::cout << " = ";
					}
					std::cout << " t = " << t << " --> the null hypothesis H0 is accepted." << std::endl;
				}
			}
		}
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::TestRegressionCoefficientsSignificanceTransformed] error testing significance of regression coefficients.");
		throw myException;
	}
}

// perform significance test
bool NuTo::MultipleLinearRegression::PerformRegressionSignificanceTest(double rTotalSumOfSquares, double rRegressionSumOfSquares, double rErrorSumOfSquares, int rNumSupportPoints, int rNumRegressorVariables, double rAlpha) const
{
	try
	{
		// calculate mean squares
		double MSr = rRegressionSumOfSquares / static_cast<double>(rNumRegressorVariables);
		double MSe = rErrorSumOfSquares / static_cast<double>(rNumSupportPoints - rNumRegressorVariables - 1);

		// calculate F-statistic
		double F0 = MSr/MSe;

		// create F-distribution
		boost::math::fisher_f FDistribution(rNumRegressorVariables, rNumSupportPoints - rNumRegressorVariables - 1);
		double F =  boost::math::quantile(boost::math::complement(FDistribution, rAlpha));

		// prepare conclusion
		bool result;
		if(F0 > F)
		{
			result = true;
		}
		else
		{
			result = false;
		}

		// prepare output
		if(this->mVerboseLevel > 0)
		{
			std::cout << std::endl;
			std::cout << "Multiple linear regression: regression significance test" << std::endl;
			std::stringstream doubleStream;
			doubleStream << std::right << std::scientific << std::setw(20) << std::setprecision(10) << std::showpos;
			std::stringstream intStream;
			intStream << std::right << std::setw(20);

			if(this->mVerboseLevel > 1)
			{
				// calculate statistics
				std::cout << "  |  Source of variation " << "  |  ";
				std::cout << "   Sum of Squares   " << "  |  ";
				std::cout << " Degrees of Freedom " << "  |  ";
				std::cout << "    Mean Square     " << "  |  ";
				std::cout << "        F0          " << "  |" << std::endl;


				std::cout << "  |    Regression        " << "  |  ";
				doubleStream << rRegressionSumOfSquares;
				std::cout << doubleStream.str() << "  |  ";
				doubleStream.str("");
				intStream << rNumRegressorVariables;
				std::cout << intStream.str() << "  |  ";
				intStream.str("");
				doubleStream << std::setw(20) << MSr;
				std::cout << doubleStream.str() << "  |  ";
				doubleStream.str("");
				doubleStream << std::setw(20) << F0;
				std::cout << doubleStream.str() << "  |  " << std::endl;
				doubleStream.str("");

				std::cout << "  |    Error or Residual " << "  |  ";
				doubleStream << std::setw(20) << rErrorSumOfSquares;
				std::cout << doubleStream.str() << "  |  ";
				doubleStream.str("");
				intStream << std::setw(20) << rNumSupportPoints - rNumRegressorVariables - 1;
				std::cout << intStream.str() << "  |  ";
				intStream.str("");
				doubleStream << std::setw(20) << MSe;
				std::cout << doubleStream.str() << "  |  ";
				doubleStream.str("");
				std::cout << "                    " << "  |  " << std::endl;

				std::cout << "  |    Total             " << "  |  ";
				doubleStream << std::setw(20) << rTotalSumOfSquares;
				std::cout << doubleStream.str() << "  |  ";
				doubleStream.str("");
				intStream << std::setw(20) << rNumSupportPoints - 1;
				std::cout << intStream.str() << "  |  ";
				intStream.str("");
				std::cout << "                    " << "  |  ";
				std::cout << "                    " << "  |  " << std::endl;
			}

			// test for significance
			std::cout << "  null hypothesis H0: all beta_i = 0.0, i=1...p" << std::endl;
			std::cout << "  alternative hypothesis H1: at least one beta_i != 0.0, i=1...p" << std::endl;
			std::cout << "  statistic: F_0 = " << F0;
			if(F0 > F)
			{
				std::cout << " > ";
			}
			else
			{
				if(F0 < F)
				{
					std::cout << " < ";
				}
				else
				{
					std::cout << " = ";
				}
			}
			std::cout << "F(" << rAlpha << "," << this->mSupportPoints.GetDimInput() << "," << this->mSupportPoints.GetNumSupportPoints() - this->mSupportPoints.GetDimInput() - 1 << ") = " << F << " --> The null hypothesis H0 is ";
			if(F0 > F)
			{
				std::cout << "rejected." << std::endl;
			}
			else
			{
				std::cout << "accepted." << std::endl;
			}
		}

		return result;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::PerformRegressionSignificanceTest] error testing significance of regression.");
		throw myException;
	}
}

// test the significance of a set of regressor variables/coefficients using original support point coordinates
bool NuTo::MultipleLinearRegression::TestGeneralRegressionSignificance(const NuTo::FullMatrix<int>& rTestCoefficients, double rAlpha) const
{
	try
	{
		return this->PerformGeneralRegressionSignificanceTest(rTestCoefficients, rAlpha, false);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::TestGeneralRegressionSignificance] error testing significance of a subset of regression coefficients.");
		throw myException;
	}
}

// test the significance of a set of regressor variables/coefficients using transformed support point coordinates
bool NuTo::MultipleLinearRegression::TestGeneralRegressionSignificanceTransformed(const NuTo::FullMatrix<int>& rTestCoefficients, double rAlpha) const
{
	try
	{
		return this->PerformGeneralRegressionSignificanceTest(rTestCoefficients, rAlpha, true);
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::TestGeneralRegressionSignificanceTransformed] error testing significance of a subset of regression coefficients.");
		throw myException;
	}
}

// test the significance of a set of regressor variables/coefficients
/*!
 * The contribution of a subset of regressor variables to the model can be investigated using the general regression significance test
 * (''extra sum of squares'' method). Let the set of \f$p\f$ regressor variables \f$\boldsymbol{x}\f$ partitioned into a subset of
 * \f$r\f$ regressor variables, \f$\boldsymbol{x}^{(1)}\f$, for which the contribution is examined and a subset of \f$p-r\f$ regressor variables,
 * \f$\boldsymbol{x}^{(2)}\f$ which are assumed to be included in the model. As a result, the multiple linear regression can be rewritten as
 * \f[
 *  \hat{y}^{(F)} = \beta_0 + \sum_{i=1}^r \beta^{(1)}_i x^{(1)}_i + \sum_{i=1}^{p-r} \beta^{(2)}_i x^{(2)}_i,
 * \f]
 * where \f$\boldsymbol{\beta^{(1)}}, \boldsymbol{\beta^{(2)}}\f$ are the corresponding subsets of regressor coefficients, and
 * \f$\hat{y}^{(F)}\f$ is the response of the full model. The error sum of squares and the regression sum of squares of the full model
 * are defined as
 * \f{align*}{
 *  SS_e &= \sum_{i=1}^n y_i^2 - \sum_{i=1}^n y_i \hat{y}^{(F)}_i\\
 *  SS_r(\beta_0, \boldsymbol{\beta^{(1)}}, \boldsymbol{\beta^{(2)})} &= \sum_{i=1}^n y_i \hat{y}^{(F)}_i- \dfrac{1}{n}\left[\sum_{i=1}^n y_i\right]^2.
 * \f}
 * In order to test the influence of \f$\boldsymbol{x}^{(1)}\f$ the following hypothesis is tested
 * \f{align*}{
 *  \text{null hypothesis }H_0: \beta^{(1)}_i &= 0.0 & i = 1 \ldots r\\
 *  \text{alternative hypothesis } H_1: \beta^{(1)}_i &\neq 0.0 & \text{for at least one }i.
 * \f}
 * Assuming that the null hypothesis \f$H_0\f$ is true, a reduced model is obtained
 * \f[
 *  \hat{y}^{(R)} = \hat{\beta}_0 + \sum_{i=1}^{p-r} \hat{\beta}^{(2)}_i x^{(2)}_i,
 * \f]
 * where \f$\hat{\beta}_0, \boldsymbol{\hat{\beta}^{(2)}}\f$ are the regression coefficients of the reduced model, which are obtaiend using the method of least squares,
 * and \f$\hat{y}^{(R)}\f$ is the response of the reduced model. The regression sum of squares of the reduced model reads
 * \f[
 *  SS_r(\beta_0, \boldsymbol{\beta^{(2)})} = \sum_{i=1}^n y_i \hat{y}^{(R)}_i- \dfrac{1}{n}\left[\sum_{i=1}^n y_i\right]^2.
 * \f]
 * The null hypothesis \f$H_0\f$ is rejected, if
 * \f[
 *  F_0 = \dfrac{SS_r(\boldsymbol{\beta^{(1)})} | \beta_0, \boldsymbol{\beta^{(2)}}) (n-p)}{SS_e r} = \dfrac{[SS_r(\beta_0, \boldsymbol{\beta^{(1)}}, \boldsymbol{\beta^{(2)})} - SS_r(\hat{\beta}_0, \boldsymbol{\hat{\beta}^{(2)})}] (n-p)}{SS_e r} > F_{\alpha,r,n-p},
 * \f]
 * where \f$SS_r(\boldsymbol{\beta^{(1)})} | \beta_0, \boldsymbol{\beta^{(2)}})\f$ is the increase in the regression sum of squares
 * due to the inclusion of \f$\boldsymbol{x^{(1)}}\f$, and \f$F_{\alpha,r,n-p}\f$ is the \f$(1-\alpha)\f$-quantile of a F-distribution with \f$(r,n-p)\f$ degrees of freedom.
 */
bool NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest(const NuTo::FullMatrix<int>& rTestCoefficients, double rAlpha, bool rTransformedFlag) const
{
	// test input data
	if(rTestCoefficients.GetNumColumns() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] number of columns in test coefficient matrix must be 1.");
	}
	if(rAlpha <= 0.0 || rAlpha >= 1.0)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] the confidence level must be larger than zero and smaller than one.");
	}
	// test model/support point compatibility
	if(this->mSupportPoints.GetDimOutput() != 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] the number of output coordinates of the support points must be equal to one.");
	}
	if(this->mSupportPoints.GetDimInput() + 1 != this->mCoefficients.GetNumRows())
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] invalid number of input coordinates of the support points. Build  model first.");
	}
	int numCoefficients = this->mCoefficients.GetNumRows();
	int numSamples = this->mSupportPoints.GetNumSupportPoints();
	if(numSamples - numCoefficients < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] invalid number of support points.");
	}

	// get test coefficients
	std::vector<bool> coefficientFlag(numCoefficients, true);
	int numTestCoefficients=0;
	for(int coefficient = 0; coefficient < rTestCoefficients.GetNumRows(); coefficient++)
	{
		int testCoefficient = rTestCoefficients(coefficient,0);
		if(testCoefficient < 1 || testCoefficient >= this->mCoefficients.GetNumRows())
		{
			throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] invalid coefficient for significance test.");
		}
		if(coefficientFlag[testCoefficient] == true)
		{
			coefficientFlag[testCoefficient] = false;
			numTestCoefficients++;
		}
	}
	if(numTestCoefficients < 1)
	{
		throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] invalid number of test coefficients.");
	}

	try
	{
		// get error sum of squares
		double SSeFull;
		if(rTransformedFlag)
		{
			SSeFull = this->GetErrorSumOfSquaresTransformed();
		}
		else
		{
			SSeFull = this->GetErrorSumOfSquares();
		}

		// get regression sum of squares (all regression coefficients)
		double SSrFull;
		if(rTransformedFlag)
		{
			SSrFull = this->GetRegressionSumOfSquaresTransformed();
		}
		else
		{
			SSrFull = this->GetRegressionSumOfSquares();
		}

		// extract reduced support point input
		const NuTo::FullMatrix<double>& fullSupportPointInput = this->mSupportPoints.GetTransformedSupportPointsInput();
		NuTo::FullMatrix<double> reducedSupportPointInput(0, this->mSupportPoints.GetNumSupportPoints());
		for(int coefficient = 1; coefficient < this->mCoefficients.GetNumRows(); coefficient++)
		{
			if(coefficientFlag[coefficient] == true)
			{
				reducedSupportPointInput.AppendRows(fullSupportPointInput.GetRow(coefficient-1));
			}
		}
		// create reduced model
		NuTo::MultipleLinearRegression reducedModel;
		reducedModel.SetSupportPoints(reducedSupportPointInput.GetNumRows(), 1, reducedSupportPointInput, this->mSupportPoints.GetTransformedSupportPointsOutput());
		reducedModel.BuildTransformation();
		reducedModel.Build();

		// calculate output of the reduced model
		NuTo::FullMatrix<double> reducedModelOutput;
		reducedModel.Solve(reducedSupportPointInput, reducedModelOutput);
		if(!rTransformedFlag)
		{
			this->mSupportPoints.TransformForwardOutput(reducedModelOutput);
		}


		// get regression sum of squares of the reduced model
		double SSrReduced = 0.0;
		if(rTransformedFlag)
		{
			// get data
			const double* reducedModelOutputData = reducedModelOutput.mEigenMatrix.data();
			const double* supportPointOutputData = this->mSupportPoints.GetTransformedSupportPointsOutput().mEigenMatrix.data();

			// calculate sums
			double sum1 = 0.0;
			double sum2 = 0.0;
			for(int sample = 0; sample < numSamples; sample++)
			{
				sum1 += supportPointOutputData[sample] * reducedModelOutputData[sample];
				sum2 += supportPointOutputData[sample];
			}

			// calculate regression sum of squares
			SSrReduced = sum1 - sum2 * sum2 / static_cast<double>(numSamples);
		}
		else
		{
			// get data
			const double* reducedModelOutputData = reducedModelOutput.mEigenMatrix.data();
			const double* supportPointOutputData = this->mSupportPoints.GetOrigSupportPointsOutput().mEigenMatrix.data();

			// calculate sums
			double sum1 = 0.0;
			double sum2 = 0.0;
			for(int sample = 0; sample < numSamples; sample++)
			{
				sum1 += supportPointOutputData[sample] * reducedModelOutputData[sample];
				sum2 += supportPointOutputData[sample];
			}

			// calculate regression sum of squares
			SSrReduced = sum1 - sum2 * sum2 / static_cast<double>(numSamples);
		}

		// calculate statistics
		double F0 = (SSrFull - SSrReduced) * (numSamples - numCoefficients) /(SSeFull * numTestCoefficients);

		// calculate quantile of F-distribution
		boost::math::fisher_f FDistribution(numTestCoefficients, numSamples - numCoefficients);
		double F =  boost::math::quantile(boost::math::complement(FDistribution, rAlpha));

		// prepare return value
		bool result;
		if(F0 > F)
		{
			result = true;
		}
		else
		{
			result = false;
		}

		// output
		if(this->mVerboseLevel > 0)
		{
			std::cout << std::endl;
			std::cout << "Multiple linear regression: general regression significance test" << std::endl;
			std::cout << "  null hypothesis H0: ";
			for(int coefficient = 0; coefficient < numCoefficients; coefficient++)
			{
				if(coefficientFlag[coefficient] == false)
				{
					std::cout << "beta_" << coefficient << " = ";
				}
			}
			std::cout << "0.0" << std::endl;
			std::cout << "  alternative hypothesis H1: ";
			bool printOrFlag = false;
			for(int coefficient = 0; coefficient < numCoefficients; coefficient++)
			{
				if(coefficientFlag[coefficient] == false)
				{
					if(printOrFlag)
					{
						std::cout << " or ";
					}
					else
					{
						printOrFlag = true;
					}
					std::cout << "beta_" << coefficient << " != 0.0";
				}
			}
			std::cout << std::endl;
			std::cout << "  statistic: F_0 = " << F0;
			if(F0 > F)
			{
				std::cout << " > ";
			}
			else
			{
				if(F0 < F)
				{
					std::cout << " < ";
				}
				else
				{
					std::cout << " = ";
				}
			}
			std::cout << "F(" << rAlpha << "," << numTestCoefficients << "," << numSamples - numCoefficients << ") = " << F << " --> The null hypothesis H0 is ";
			if(result)
			{
				std::cout << "rejected." << std::endl;
			}
			else
			{
				std::cout << "accepted." << std::endl;
			}
		}

		// return
		return result;
	}
	catch(NuTo::Exception& e)
	{
		NuTo::MetamodelException myException(e.ErrorMessage());
		myException.AddMessage("[NuTo::MultipleLinearRegression::PerformGeneralRegressionSignificanceTest] error performing significance test on a subset of regression coefficients.");
		throw myException;
	}
}


#ifdef ENABLE_SERIALIZATION
// restore the object from a file
void NuTo::MultipleLinearRegression::Restore (const std::string &filename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

        // open file
        std::ifstream ifs (filename.c_str(), std::ios_base::binary );
        if (!ifs.is_open())
        {
            throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Restore] error opening file");
        }

        // read date
        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if (typeIdString != this->GetTypeId())
            {
                throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if (typeIdString != this->GetTypeId())
            {
                throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if (typeIdString != this->GetTypeId())
            {
                throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Restore]Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw NuTo::MetamodelException( "[NuTo::MultipleLinearRegression::Restore]File type not implemented" );
        }

        // close file
        ifs.close();
    }
    catch (NuTo::MetamodelException &e)
    {
        throw e;
    }
    catch (std::exception &e)
    {
        throw NuTo::MetamodelException(e.what());
    }
    catch ( ... )
    {
        throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Restore] Unhandled exception.");
    }
}

// save the object to a file
void NuTo::MultipleLinearRegression::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if (!ofs.is_open())
        {
            throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Save] error opening file.");
        }

        // write data
        std::string typeIdString ( this->GetTypeId() );
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Save] File type not implemented.");
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::MultipleLinearRegression::Save]File save exception in boost - " ) +std::string ( e.what() ) );
        throw NuTo::MetamodelException(s);
    }
    catch ( MathException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw NuTo::MetamodelException(e.what());
    }
    catch ( ... )
    {
        throw NuTo::MetamodelException("[NuTo::MultipleLinearRegression::Save] Unhandled exception.");
    }
}

// serializes the class
template void NuTo::MultipleLinearRegression::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::MultipleLinearRegression::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::MultipleLinearRegression::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::MultipleLinearRegression::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::MultipleLinearRegression::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::MultipleLinearRegression::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::MultipleLinearRegression::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize MultipleLinearRegression" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Metamodel)
       & BOOST_SERIALIZATION_NVP(mCoefficients);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize MultipleLinearRegression" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::MultipleLinearRegression)
#endif  // ENABLE_SERIALIZATION
