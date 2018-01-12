#pragma once

#include <functional>
#include "mechanics/elements/ContinuumBoundaryElement.h"
#include "mechanics/IGA/NURBSCurve.h"


namespace NuTo
{
struct ContactGap
{
    //! @brief function defining the gap prior to contact
	std::function<Eigen::VectorXd(Eigen::VectorXd)> GapFunction;

	//! @brief and its first and second derivative
    std::function<Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>(Eigen::VectorXd)> GapFunctionDerivative1;
    std::function<Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>(Eigen::VectorXd)> GapFunctionDerivative2;

    //! @brief penalty parameter for the constitutive contact law
    double PenaltyParameter;
};

struct ContactGapNURBS
{
	//! @brief coordinates of the points on the master surface prior to penetration
	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> CoordinatesAndIDs;

	//! @brief degree of the polynomial
	int Degree;

	//! @brief penalty parameter for the constitutive contact law
	double PenaltyParameter;
};
//NuTo::ContactGap DefineContactGap(std::function<Eigen::VectorXd(Eigen::VectorXd)> rGapFunction,
//	    std::function<Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>(Eigen::VectorXd)> rGapFunctionDerivative1,
//	    std::function<Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic>(Eigen::VectorXd)> rGapFunctionDerivative2,
//		double* rPenaltyParameter = nullptr)
//{
//	if (rPenaltyParameter) {
//		return NuTo::ContactGap({rGapFunction, rGapFunctionDerivative1, rGapFunctionDerivative2, *rPenaltyParameter});
//	} else {
//		return NuTo::ContactGap({rGapFunction, rGapFunctionDerivative1, rGapFunctionDerivative2, 1.});
//	}
//};

template <int TDim>
class ContinuumBoundaryElementContact : public ContinuumBoundaryElement<TDim>
{
public:
    ContinuumBoundaryElementContact(const ContinuumElement<TDim>& rBaseElement,
            				const IntegrationTypeBase& integrationType,
							int rSurfaceId);

    virtual ~ContinuumBoundaryElementContact() = default;

    //! @brief getter for the ContactGap
    void SetGeometryRigidBody(ContactGap& rGapProblem)
    {
    	mGapProblem.GapFunction = rGapProblem.GapFunction;
    	mGapProblem.GapFunctionDerivative1 = rGapProblem.GapFunctionDerivative1;
    	mGapProblem.GapFunctionDerivative2 = rGapProblem.GapFunctionDerivative2;
    	mGapProblem.PenaltyParameter = rGapProblem.PenaltyParameter;

    	mPenaltyParameter = rGapProblem.PenaltyParameter;
    }

    //! @brief getter for the ContactGapNURBS
    void SetGeometryNURBS(ContactGapNURBS& rGapProblemNURBS)
    {
    	mGapProblemNURBS.CoordinatesAndIDs = rGapProblemNURBS.CoordinatesAndIDs;
    	mGapProblemNURBS.Degree = rGapProblemNURBS.Degree;
    	mGapProblemNURBS.PenaltyParameter = rGapProblemNURBS.PenaltyParameter;
    }

    void SetNURBScurve(NURBSCurve rNURBScurve)
    {
    	mNURBScurve = rNURBScurve;

    	mNURBSinterpolation = true;

    	mPenaltyParameter = mGapProblemNURBS.PenaltyParameter;
    }

//    double GetPenaltyParameter() const
//    {
//        return mPenaltyParameter;
//    }

    //! @brief getter for the penalty parameter
//    double GetPenaltyParameter() const
//    {
//        return mPenaltyParameter;
//    }

    //! @brief setter for the penalty parameter
//    void SetPenaltyParameter(double rPenaltyParameter)
//    {
//        mPenaltyParameter = rPenaltyParameter;
//    }

protected:
    void CalculateElementOutputInternalGradient(BlockFullVector<double>& rInternalGradient,
                                                EvaluateDataContinuumBoundary<TDim>& rData,
                                                const ConstitutiveInputMap& constitutiveInput,
                                                const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const override;

    void CalculateElementOutputHessian0(BlockFullMatrix<double>& rHessian0, EvaluateDataContinuumBoundary<TDim>& rData,
                                        const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const;

    void CalculateElementOutputHessianAxSy0(BlockFullMatrix<double>& rHessian0, EvaluateDataContinuumBoundary<TDim>& rData,
                                        const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const
    {
    	throw NuTo::Exception(__PRETTY_FUNCTION__, "Only implemented for 2D AXISYMMETRIC.");
    }

    void CalculateNmatrixNormal(EvaluateDataContinuumBoundary<TDim>& rData, Eigen::MatrixXd& rNmatrixNormal, double& rGap, int rTheIP) const;

    std::pair<Eigen::VectorXd, double> CalculateNormalToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const;

    Eigen::VectorXd ProjectionToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const;

    Eigen::VectorXd GetGlobalIntegrationPointCoordinatesAfterDeformation(EvaluateDataContinuumBoundary<TDim>& rData, int rIpNum) const;

    void InterpolateNURBScurce();

    ContactGap mGapProblem;

    ContactGapNURBS mGapProblemNURBS;

    NURBSCurve mNURBScurve;

    double mPenaltyParameter;

    bool mNURBSinterpolation;

};

} /* namespace NuTo */
