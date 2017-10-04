/*
 * ContinuumBoundaryElementBase.cpp
 *
 *  Created on: 5 Mar 2016
 *      Author: vhirtham
 */
#include <iostream>

#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/dofSubMatrixStorage/BlockFullMatrix.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/elements/ContinuumBoundaryElementContact.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/elements/ElementOutputBase.h"
#include "mechanics/elements/ElementOutputIpData.h"
#include "mechanics/elements/EvaluateDataContinuumBoundary.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/interpolationtypes/InterpolationBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/sections/Section.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

using namespace NuTo;

template <int TDim>
NuTo::ContinuumBoundaryElementContact<TDim>::ContinuumBoundaryElementContact(const ContinuumElement<TDim>& rBaseElement,
                                                               const IntegrationTypeBase& integrationType,
                                                               int rSurfaceId)
    : ContinuumBoundaryElement<TDim>(rBaseElement, integrationType, rSurfaceId)
{

}

template <int TDim>
void NuTo::ContinuumBoundaryElementContact<TDim>::CalculateNmatrixNormal(
        EvaluateDataContinuumBoundary<TDim>& rData, Eigen::MatrixXd& rNmatrixNormal, double& rGap, int rTheIP) const
{
	std::cout << "Bin here: CalculateNmatrixNormal" << std::endl;
	Eigen::VectorXd GlobalCoordinateTheIP;
	Eigen::Vector3d GlobalCoordinateTheIP3D;

	rNmatrixNormal.setZero(TDim, rData.mN.at(Node::eDof::DISPLACEMENTS).cols());

//	Eigen::MatrixXd rDataN(rData.mN.at(Node::eDof::COORDINATES));
//	std::cout << "--- CalculateNmatrixNormal: rData" << rDataN.eval() << std::endl;


	// GlobalCoordinateTheIP = rData.mN.at(Node::eDof::COORDINATES) * rData.mNodalValues.at(Node::eDof::COORDINATES);
	//==> check above calculation of the globals slave coordinate of rTheIp
	// no, it is wrong cause rData does not contains mN for eDof::COORDINATES
	GlobalCoordinateTheIP3D = this->ContinuumBoundaryElement<TDim>::GetGlobalIntegrationPointCoordinates(rTheIP);
	GlobalCoordinateTheIP = GlobalCoordinateTheIP3D.segment(0, TDim);

	std::cout << "--- CalculateNmatrixNormal: " << GlobalCoordinateTheIP << std::endl;

	std::pair<Eigen::VectorXd, double> NormalGap = CalculateNormalToRigidBody(GlobalCoordinateTheIP);
	Eigen::VectorXd normal = NormalGap.first;

	std::cout << "--- CalculateNmatrixNormal: get normal =" << normal.transpose() << std::endl;

	assert(rData.mN.at(Node::eDof::DISPLACEMENTS).rows() == normal.rows());

//	for (int i = 0; i < normal.rows(); ++i) {
//		std::cout << "--- CalculateNmatrixNormal: in the for loop, i = "<< i << " normal.rows() = " << normal.rows() << std::endl;
//		std::cout << "--- : normal(i) = " << normal(i) << std::endl;
//		std::cout << "--- : rData.mN.at(Node::eDof::DISPLACEMENTS).row(i) = " << rData.mN.at(Node::eDof::DISPLACEMENTS).row(i) << std::endl;
//
//		rNmatrixNormal.row(i) = normal(i) * rData.mN.at(Node::eDof::DISPLACEMENTS).row(i);
//	}

	Eigen::MatrixXd matrixN, NmatrixNormalT;

	matrixN = rData.mN.at(Node::eDof::DISPLACEMENTS);
	NmatrixNormalT = matrixN.transpose() * normal;

	rNmatrixNormal = NmatrixNormalT.transpose();

	std::cout << "--- CalculateNmatrixNormal: NmatrixNormal " << std::endl;
	std::cout << rNmatrixNormal << std::endl;

	rGap = NormalGap.second;

	std::cout << "--- CalculateNmatrixNormal: rGap = "<< rGap << std::endl;
}

//template <int TDim>
//void NuTo::ContinuumBoundaryElementContact<TDim>::CalculateElementOutputs(
//        std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>>& rElementOutput,
//        EvaluateDataContinuumBoundary<TDim>& rData, int rTheIP, const ConstitutiveInputMap& constitutiveInput,
//        const ConstitutiveOutputMap& constitutiveOutput) const
//{
//	std::cout << "Bin hier: CalculateElementOutputs" << std::endl;

//this->ContinuumBoundaryElement<TDim>::CalculateElementOutputs(rElementOutput, rData, rTheIP, constitutiveInput, constitutiveOutput);

//    rData.mDetJxWeightIPxSection =
//            CalculateDetJxWeightIPxSection(rData.mDetJacobian, rTheIP); // formerly known as "factor"
//
//    for (auto it : rElementOutput)
//    {
//        switch (it.first)
//        {
//        case Element::eOutput::INTERNAL_GRADIENT:
//            UpdateAlphaGradientDamage(rData, constitutiveInput, constitutiveOutput);
//            CalculateElementOutputInternalGradient(it.second->GetBlockFullVectorDouble(), rData, constitutiveInput,
//                                                   constitutiveOutput, rTheIP);
//            break;
//
//        case Element::eOutput::HESSIAN_0_TIME_DERIVATIVE:
//            UpdateAlphaGradientDamage(rData, constitutiveInput, constitutiveOutput);
//            CalculateElementOutputHessian0(it.second->GetBlockFullMatrixDouble(), rData, constitutiveOutput, rTheIP);
//            break;
//
//        case Element::eOutput::HESSIAN_1_TIME_DERIVATIVE:
//            break;
//
//        case Element::eOutput::HESSIAN_2_TIME_DERIVATIVE:
//            break;
//
//        case Element::eOutput::LUMPED_HESSIAN_2_TIME_DERIVATIVE:
//            break;
//
//        case Element::eOutput::UPDATE_STATIC_DATA:
//        case Element::eOutput::UPDATE_TMP_STATIC_DATA:
//            break;
//        case Element::eOutput::IP_DATA:
//            CalculateElementOutputIpData(it.second->GetIpData(), constitutiveOutput, rTheIP);
//            break;
//        case Element::eOutput::GLOBAL_ROW_DOF:
//        case Element::eOutput::GLOBAL_COLUMN_DOF:
//            break;
//        default:
//            throw Exception(__PRETTY_FUNCTION__, "element output not implemented.");
//        }
//    }
//}

template <int TDim>
void NuTo::ContinuumBoundaryElementContact<TDim>::CalculateElementOutputInternalGradient(
        BlockFullVector<double>& rInternalGradient, EvaluateDataContinuumBoundary<TDim>& rData,
        const ConstitutiveInputMap& constitutiveInput, const ConstitutiveOutputMap& constitutiveOutput,
        int rTheIP) const
{
    for (auto dofRow : ContinuumBoundaryElement<TDim>::mInterpolationType->GetActiveDofs())
    {
        switch (dofRow)
        {
        case Node::eDof::DISPLACEMENTS:
        {
        	std::cout << "Bin here: CalculateElementOutputInternalGradient" << std::endl;
        	Eigen::MatrixXd NmatrixNormal;
        	double gap;

        	CalculateNmatrixNormal(rData, NmatrixNormal, gap, rTheIP);

        	//==> update rInternalGradient to the contact force
        	if (gap < 0.) {
        		std::cout << "*** *** CONTRIBUTE TO INTERNAL GRADIENT *** ***" << std::endl;
            	rInternalGradient[dofRow] +=
            			rData.mDetJxWeightIPxSection * NmatrixNormal.transpose() * mGapProblem.PenaltyParameter * gap;
			}
        	break;
        }

        case Node::eDof::NONLOCALEQSTRAIN:
        {
            const auto& localEqStrain = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::LOCAL_EQ_STRAIN).get());
            const auto& nonlocalEqStrain = *static_cast<ConstitutiveScalar*>(
                    constitutiveInput.at(Constitutive::eInput::NONLOCAL_EQ_STRAIN).get());

            rInternalGradient[dofRow] += rData.mDetJxWeightIPxSection * rData.m1DivAlpha *
                                         rData.mN.at(Node::eDof::NONLOCALEQSTRAIN).transpose() *
                                         (nonlocalEqStrain[0] - localEqStrain[0]);
            break;
        }

        case Node::eDof::RELATIVEHUMIDITY:
        {
            const auto& internalGradientRH_Boundary_N = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::INTERNAL_GRADIENT_RELATIVE_HUMIDITY_BOUNDARY_N).get());
            rInternalGradient[dofRow] +=
                    rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * internalGradientRH_Boundary_N;
            break;
        }
        case Node::eDof::WATERVOLUMEFRACTION:
        {
            const auto& internalGradientWV_Boundary_N = *static_cast<ConstitutiveScalar*>(
                    constitutiveOutput.at(Constitutive::eOutput::INTERNAL_GRADIENT_WATER_VOLUME_FRACTION_BOUNDARY_N)
                            .get());
            rInternalGradient[dofRow] +=
                    rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() * internalGradientWV_Boundary_N;
            break;
        }
        default:
            throw Exception(__PRETTY_FUNCTION__, "Element output INTERNAL_GRADIENT for " +
                                                                  Node::DofToString(dofRow) + " not implemented.");
        }
    }
}


template <int TDim>
void NuTo::ContinuumBoundaryElementContact<TDim>::CalculateElementOutputHessian0(
        BlockFullMatrix<double>& rHessian0, EvaluateDataContinuumBoundary<TDim>& rData,
        const ConstitutiveOutputMap& constitutiveOutput, int rTheIP) const
{
    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    for (auto dofRow : ContinuumBoundaryElement<TDim>::mInterpolationType->GetActiveDofs())
    {
        for (auto dofCol : ContinuumBoundaryElement<TDim>::mInterpolationType->GetActiveDofs())
        {
            if (!ContinuumBoundaryElement<TDim>::GetConstitutiveLaw(rTheIP).CheckDofCombinationComputable(dofRow, dofCol, 0))
                continue;
            auto& hessian0 = rHessian0(dofRow, dofCol);
            switch (Node::CombineDofs(dofRow, dofCol))
            {

            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
			{
            	std::cout << "Bin here: CalculateElementOutputHESSIAN0" << std::endl;
            	Eigen::MatrixXd NmatrixNormal;
            	double gap;

            	CalculateNmatrixNormal(rData, NmatrixNormal, gap, rTheIP);

            	//==> update stiffness due to the contact force
            	if (gap < 0.) {
            		std::cout << "*** *** CONTRIBUTE TO HASSIAN *** ***" << std::endl;
            		hessian0 += rData.mDetJxWeightIPxSection * NmatrixNormal.transpose() *
            				mGapProblem.PenaltyParameter * NmatrixNormal;
            	}
            	break;
			}
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::NONLOCALEQSTRAIN):
                break;

            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::DISPLACEMENTS):
            {
                const auto& tangentLocalEqStrainStrain = *static_cast<ConstitutiveVector<VoigtDim>*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN).get());
                hessian0 -= rData.mDetJxWeightIPxSection * rData.m1DivAlpha * rData.mN.at(dofRow).transpose() *
                            tangentLocalEqStrainStrain.transpose() * rData.mB.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::NONLOCALEQSTRAIN, Node::eDof::NONLOCALEQSTRAIN):
            {
                hessian0 += rData.mN.at(dofRow).transpose() * rData.mN.at(dofRow) * rData.mDetJxWeightIPxSection *
                            rData.m1DivAlpha;
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::RELATIVEHUMIDITY):
            {
                const auto& internalGradientRH_dRH_Boundary_NN_H0 = *static_cast<ConstitutiveScalar*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_INTERNAL_GRADIENT_RH_D_RH_BOUNDARY_NN_H0).get());
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() *
                            internalGradientRH_dRH_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::WATERVOLUMEFRACTION):
            {
                const auto& internalGradientWV_dWV_Boundary_NN_H0 = *static_cast<ConstitutiveScalar*>(
                        constitutiveOutput.at(Constitutive::eOutput::D_INTERNAL_GRADIENT_WV_D_WV_BOUNDARY_NN_H0).get());
                hessian0 += rData.mDetJxWeightIPxSection * rData.mN.at(dofRow).transpose() *
                            internalGradientWV_dWV_Boundary_NN_H0 * rData.mN.at(dofCol);
                break;
            }
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::RELATIVEHUMIDITY):
                break;
            default:
            /*******************************************************\
            |         NECESSARY BUT UNUSED DOF COMBINATIONS         |
            \*******************************************************/
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::RELATIVEHUMIDITY):
            case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::WATERVOLUMEFRACTION):
            case Node::CombineDofs(Node::eDof::RELATIVEHUMIDITY, Node::eDof::DISPLACEMENTS):
            case Node::CombineDofs(Node::eDof::WATERVOLUMEFRACTION, Node::eDof::DISPLACEMENTS):
                continue;
                throw Exception(__PRETTY_FUNCTION__, "Element output HESSIAN_0_TIME_DERIVATIVE for "
                                                              "(" + Node::DofToString(dofRow) +
                                                                      "," + Node::DofToString(dofCol) +
                                                                      ") not implemented.");
            }
        }
    }
}

template <>
std::pair<Eigen::VectorXd, double> NuTo::ContinuumBoundaryElementContact<1>::CalculateNormalToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const
{
	throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
		                                     "] normal to a 1D rigid body!.");
}

template <>
std::pair<Eigen::VectorXd, double> NuTo::ContinuumBoundaryElementContact<2>::CalculateNormalToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const
{
	std::cout << "Bin here: CalculateNormalToRigidBody" << std::endl;
	// calculate the projection of rSlaveCoordinate to the rigid surface
	Eigen::VectorXd ProjectionToMaster = ProjectionToRigidBody(rSlaveCoordinates);

	std::cout << "CalculateNormalToRigidBody: projection to master = " << ProjectionToMaster.transpose() << std::endl;

	Eigen::VectorXd tangent = mGapProblem.GapFunctionDerivative1(ProjectionToMaster)(0,0);

	std::cout << "- - -: tangent = " << tangent.transpose() << std::endl;

	if(tangent.rows() != 2)
		throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
						                                     "] the normal to a 2D rigid body is expected.");

	Eigen::Vector2d normal;
    normal(0) = tangent(1);
    normal(1) = -tangent(0);

    normal.normalize();

    Eigen::VectorXd projectionVector(rSlaveCoordinates - ProjectionToMaster);
    Eigen::VectorXd normalXd;

	normalXd = normal.segment(0, 2);

    return std::make_pair(normalXd, projectionVector.dot(normal));
}

template <>
std::pair<Eigen::VectorXd, double> NuTo::ContinuumBoundaryElementContact<3>::CalculateNormalToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const
{
	throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
		                                     "] normal to a 3D rigid body is not implemented!.");
}

template <int TDim>
Eigen::VectorXd NuTo::ContinuumBoundaryElementContact<TDim>::ProjectionToRigidBody(Eigen::VectorXd &rSlaveCoordinates) const
{
	std::cout << "Bin here: ProjectionToRgidBody" << std::endl;
	double tol = 1.e-8;
	double error = 1.;
	int maxNumIter = 20;
	int numIter = 0;

	int numPrimes = TDim - 1;

	// define the residual and the parameter vector (the latter is unknown which parameterizes the rigid surface)
	Eigen::VectorXd dprime(numPrimes), parameters(numPrimes), coordinatesMaster;

	// initialize parameters
	parameters.setZero(numPrimes);

	try {
		coordinatesMaster = mGapProblem.GapFunction(parameters);
	} catch (...) {
		throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
				                                     "] during initialization: the rigid surface is undefined at zero parameters!.");
	}

	if(mGapProblem.GapFunctionDerivative1(parameters).cols() != numPrimes)
		throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
								                                     "] the derivative1 matrix should have" + std::to_string(numPrimes) + " columns.");
	if(mGapProblem.GapFunctionDerivative2(parameters).cols() != numPrimes || mGapProblem.GapFunctionDerivative2(parameters).rows() != numPrimes)
		throw Exception(std::string("[") + __PRETTY_FUNCTION__ +
								                                     "] the derivative2 matrix should have" + std::to_string(numPrimes) + " columns and rows.");

	Eigen::VectorXd projectionVector;

	std::cout << "- - - : prior to NR" << std::endl;

	// newton-raphson
	while(error > tol && numIter < maxNumIter)
	{
		// update coordinatesMaster
		coordinatesMaster = mGapProblem.GapFunction(parameters);
		projectionVector = rSlaveCoordinates - coordinatesMaster;

		dprime.setZero(numPrimes);

		Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> prime = mGapProblem.GapFunctionDerivative1(parameters);
		Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> primeprime = mGapProblem.GapFunctionDerivative2(parameters);

		// calculate residual dprime
        for(int j = 0; j < prime.cols(); j++)
        	dprime(j) += projectionVector.dot(prime(0,j));

        // calculate derivative dprimeprime
        Eigen::MatrixXd dprimeprime(numPrimes, numPrimes);
        dprimeprime.setZero(numPrimes, numPrimes);

        for(int i = 0; i < prime.cols() ; i++)
        {
            for(int j = 0; j < prime.cols() ; j++)
            {
                dprimeprime(i,j) = projectionVector.dot(primeprime(i,j)) - prime(0,i).dot(prime(0,j));
            }
        }

        // solve for increment
        Eigen::VectorXd increment = dprimeprime.colPivHouseholderQr().solve(-dprime);

        parameters += increment;

        error = dprime.norm();
        numIter++;

		std::cout << "- - - : numIter " << numIter << std::endl;
		std::cout << "- - - : error " << error << std::endl;
		std::cout << "- - - : coordMaster " << coordinatesMaster.transpose() << std::endl;
	}

	 if(numIter >= maxNumIter)
		 std::cout << "NuTo::ContinuumBoundaryElementContact::ProjectionToRigidBody]: Maximum number of Newton iterations exceeded!" << std::endl;

	std::cout << "- - - : out of ProjectionToRigidBody " << coordinatesMaster.transpose() << std::endl;

	return coordinatesMaster;
}

template class NuTo::ContinuumBoundaryElementContact<1>;
template class NuTo::ContinuumBoundaryElementContact<2>;
template class NuTo::ContinuumBoundaryElementContact<3>;
