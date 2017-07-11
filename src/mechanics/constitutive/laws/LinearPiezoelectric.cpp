#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"

#include "mechanics/constitutive/laws/LinearPiezoelectric.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "base/Logger.h"
#include "base/Exception.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

#include <boost/current_function.hpp>

NuTo::LinearPiezoelectric::LinearPiezoelectric()
    : ConstitutiveBase()
{
    mStiffness << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0;
    mRho = 0.;

    mPermittivity << 0, 0, 0, 0, 0, 0, 0, 0, 0;

    mPiezo << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;

    SetParametersValid();
}

//#ifdef ENABLE_SERIALIZATION
////! @brief serializes the class
////! @param ar         archive
////! @param version    version
// template<class Archive>
// void NuTo::LinearPiezoelectric::serialize(Archive & ar, const unsigned int version)
//{
//#ifdef DEBUG_SERIALIZATION
//    std::cout << "start serialize LinearPiezoelectric" << std::endl;
//#endif
//    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
//    & BOOST_SERIALIZATION_NVP(mE)
//    & BOOST_SERIALIZATION_NVP(mNu)
//    & BOOST_SERIALIZATION_NVP(mRho)
//    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
//#ifdef DEBUG_SERIALIZATION
//    std::cout << "finish serialize LinearPiezoelectric" << std::endl;
//#endif
//}
// BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LinearPiezoelectric)
//#endif // ENABLE_SERIALIZATION

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::LinearPiezoelectric::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<LinearPiezoelectric>>(*this);
}


NuTo::ConstitutiveInputMap
NuTo::LinearPiezoelectric::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput,
                                                 const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        case NuTo::Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            constitutiveInputMap[Constitutive::eInput::ELECTRIC_FIELD];
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            break;
        case NuTo::Constitutive::eOutput::ELECTRIC_FIELD:
            constitutiveInputMap[Constitutive::eInput::ELECTRIC_FIELD];
            break;
        // no inputs needed for:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
            //            ProcessUnhandledOutput(__PRETTY_FUNCTION__,itOutput.first);
            //            throw Exception(std::string("[")+__PRETTY_FUNCTION__+"] output object " +
            //            Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive
            //            law.");
        }
    }

    return constitutiveInputMap;
}

namespace NuTo // template specialization in same namespace as definition
{
template <>
void NuTo::LinearPiezoelectric::Evaluate<1>(const ConstitutiveInputMap& rConstitutiveInput,
                                            const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

            ConstitutiveIOBase& electricField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            electricField.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress[0] =
                    this->mStiffness(0, 0) * engineeringStrain[0] - this->mPiezo(0, 0) * electricField[0];
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

            ConstitutiveIOBase& electricField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            electricField.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress[0] =
                    this->mStiffness(0, 0) * engineeringStrain[0] - this->mPiezo(0, 0) * electricField[0];
            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_FIELD:
        {
            *itOutput.second = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            ConstitutiveIOBase& electricDisplacement = *itOutput.second;
            electricDisplacement.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

            ConstitutiveIOBase& electricField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            electricField.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            electricDisplacement[0] =
                    this->mPiezo(0, 0) * engineeringStrain[0] + this->mPermittivity(0, 0) * electricField[0];
            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = this->mStiffness(0, 0);
            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = -this->mPiezo(0, 0);
            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = this->mPiezo(0, 0);
            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = this->mPermittivity(0, 0);
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D().As3D(0);
            break;
        }
        case NuTo::Constitutive::eOutput::EXTRAPOLATION_ERROR:
            break;
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            // nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}


template <>
void NuTo::LinearPiezoelectric::Evaluate<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                            const ConstitutiveOutputMap& rConstitutiveOutput)
{
    const auto& planeState =
            *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

            ConstitutiveIOBase& eField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            eField.AssertIsVector<2>(itOutput.first, __PRETTY_FUNCTION__);

            Eigen::Matrix<double, 3, 1> electricField;
            electricField[0] = eField[0];
            electricField[1] = eField[1];
            electricField[2] = 0.;

            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStress(6);
            Eigen::VectorXd engStrain(6);
            engStrain.setZero();
            engStress.setZero();
            engStrain(0) = engineeringStrain[0];
            engStrain(1) = engineeringStrain[1];
            engStrain(5) = engineeringStrain[2];

            engStress = this->mStiffness * engStrain - (this->mPiezo).transpose() * electricField;

            engineeringStress[0] = engStress(0);
            engineeringStress[1] = engStress(1);
            engineeringStress[2] = engStress(5);

            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE: // for visualization
        {
            itOutput.second->AsEngineeringStrain3D() = rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)
                                                               ->AsEngineeringStrain2D()
                                                               .As3D(0, planeState.GetPlaneState());
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

            ConstitutiveIOBase& eField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            eField.AssertIsVector<2>(itOutput.first, __PRETTY_FUNCTION__);

            Eigen::Matrix<double, 3, 1> electricField;
            electricField[0] = eField[0];
            electricField[1] = eField[1];
            electricField[2] = 0.;


            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStress(6);
            Eigen::VectorXd engStrain(6);
            engStrain.setZero();
            engStress.setZero();
            engStrain(0) = engineeringStrain[0];
            engStrain(1) = engineeringStrain[1];
            engStrain(5) = engineeringStrain[2];

            engStress = this->mStiffness * engStrain - (this->mPiezo).transpose() * electricField;

            engineeringStress[0] = engStress(0);
            engineeringStress[1] = engStress(1);
            engineeringStress[2] = engStress(2);
            engineeringStress[3] = engStress(3);
            engineeringStress[4] = engStress(4);
            engineeringStress[5] = engStress(5);


            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_FIELD:
        {
            *itOutput.second = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            ConstitutiveIOBase& electricDisplacement = *itOutput.second;
            electricDisplacement.AssertIsVector<2>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

            ConstitutiveIOBase& eField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            eField.AssertIsVector<2>(itOutput.first, __PRETTY_FUNCTION__);

            Eigen::Matrix<double, 3, 1> electricField;
            electricField[0] = eField[0];
            electricField[1] = eField[1];
            electricField[2] = 0.;


            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStrain(6);
            engStrain.setZero();
            engStrain(0) = engineeringStrain[0];
            engStrain(1) = engineeringStrain[1];
            engStrain(5) = engineeringStrain[2];

            Eigen::VectorXd elDispl(3);

            elDispl = this->mPiezo * engStrain + this->mPermittivity * electricField;

            for (int ii = 0; ii < 2; ii++)
            {
                electricDisplacement[ii] = elDispl(ii);
            }

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = this->mStiffness(0, 0);
            tangent(0, 1) = this->mStiffness(0, 1);
            tangent(0, 2) = this->mStiffness(0, 5);
            tangent(1, 0) = this->mStiffness(1, 0);
            tangent(1, 1) = this->mStiffness(1, 1);
            tangent(1, 2) = this->mStiffness(1, 5);
            tangent(2, 0) = this->mStiffness(5, 0);
            tangent(2, 1) = this->mStiffness(5, 1);
            tangent(2, 2) = this->mStiffness(5, 5);

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 2>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = -((this->mPiezo).transpose())(0, 0);
            tangent(0, 1) = -((this->mPiezo).transpose())(0, 1);
            tangent(1, 0) = -((this->mPiezo).transpose())(1, 0);
            tangent(1, 1) = -((this->mPiezo).transpose())(1, 1);
            tangent(2, 0) = -((this->mPiezo).transpose())(5, 0);
            tangent(2, 1) = -((this->mPiezo).transpose())(5, 1);

            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<2, 3>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0, 0) = this->mPiezo(0, 0);
            tangent(0, 1) = this->mPiezo(0, 1);
            tangent(0, 2) = this->mPiezo(0, 5);
            tangent(1, 0) = this->mPiezo(1, 0);
            tangent(1, 1) = this->mPiezo(1, 1);
            tangent(1, 2) = this->mPiezo(1, 5);

            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<2, 2>(itOutput.first, __PRETTY_FUNCTION__);

            for (int ii = 0; ii < 2; ii++)
            {
                for (int jj = 0; jj < 2; jj++)
                {
                    tangent(ii, jj) = this->mPermittivity(ii, jj);
                }
            }

            break;
        }
        case NuTo::Constitutive::eOutput::EXTRAPOLATION_ERROR:
            break;
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            // nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}


template <>
void NuTo::LinearPiezoelectric::Evaluate<3>(const ConstitutiveInputMap& rConstitutiveInput,
                                            const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();

            ConstitutiveIOBase& eField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            eField.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            Eigen::Matrix<double, 3, 1> electricField = eField.CopyToEigenMatrix();

            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStress(6);
            Eigen::VectorXd engStrain(6);
            for (int ii = 0; ii < 6; ii++)
            {
                engStrain(ii) = engineeringStrain[ii];
            }

            engStress = this->mStiffness * engStrain - (this->mPiezo).transpose() * electricField;

            for (int ii = 0; ii < 6; ii++)
            {
                engineeringStress[ii] = engStress(ii);
            }

            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_FIELD:
        {
            *itOutput.second = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            ConstitutiveIOBase& electricDisplacement = *itOutput.second;
            electricDisplacement.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();

            ConstitutiveIOBase& eField = *rConstitutiveInput.at(Constitutive::eInput::ELECTRIC_FIELD);
            eField.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);
            Eigen::Matrix<double, 3, 1> electricField = eField.CopyToEigenMatrix();

            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStrain(6);
            Eigen::VectorXd elDispl(3);
            for (int ii = 0; ii < 6; ii++)
            {
                engStrain(ii) = engineeringStrain[ii];
            }

            elDispl = this->mPiezo * engStrain + this->mPermittivity * electricField;

            for (int ii = 0; ii < 3; ii++)
            {
                electricDisplacement[ii] = elDispl(ii);
            }

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<6, 6>(itOutput.first, __PRETTY_FUNCTION__);

            for (int ii = 0; ii < 6; ii++)
            {
                for (int jj = 0; jj < 6; jj++)
                {
                    tangent(ii, jj) = this->mStiffness(ii, jj);
                }
            }

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<6, 3>(itOutput.first, __PRETTY_FUNCTION__);

            for (int ii = 0; ii < 6; ii++)
            {
                for (int jj = 0; jj < 3; jj++)
                {
                    tangent(ii, jj) = -((this->mPiezo).transpose())(ii, jj);
                }
            }

            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 6>(itOutput.first, __PRETTY_FUNCTION__);

            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = 0; jj < 6; jj++)
                {
                    tangent(ii, jj) = this->mPiezo(ii, jj);
                }
            }

            break;
        }
        case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

            for (int ii = 0; ii < 3; ii++)
            {
                for (int jj = 0; jj < 3; jj++)
                {
                    tangent(ii, jj) = this->mPermittivity(ii, jj);
                }
            }

            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();
            break;
        }
        case NuTo::Constitutive::eOutput::EXTRAPOLATION_ERROR:
            break;
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            // nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

} // namespace NuTo


bool NuTo::LinearPiezoelectric::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow, NuTo::Node::eDof rDofCol,
                                                              int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
    case 2:
    {
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::ELECTRICPOTENTIAL, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::ELECTRICPOTENTIAL):
        case Node::CombineDofs(Node::eDof::ELECTRICPOTENTIAL, Node::eDof::ELECTRICPOTENTIAL):
            return true;
        default:
            return false;
        }
    }
    break;
    default:
        return false;
    }
}


bool NuTo::LinearPiezoelectric::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
    case Constitutive::eConstitutiveParameter::STIFFNESS:
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    case Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}


double NuTo::LinearPiezoelectric::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}


void NuTo::LinearPiezoelectric::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier,
                                                   double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        this->mRho = rValue;
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

Eigen::MatrixXd
NuTo::LinearPiezoelectric::GetParameterMatrixDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::STIFFNESS:
    {
        return mStiffness;
    }
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        return mPermittivity;
    }
    case Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR:
    {
        return mPiezo;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::LinearPiezoelectric::SetParameterMatrixDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier,
                                                         Eigen::MatrixXd rValue)
{
    // ConstitutiveBase::CheckParameterFullVector(rIdentifier, rValue);
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::STIFFNESS:
    {
        mStiffness = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        mPermittivity = rValue;
        break;
    }
    case Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR:
    {
        mPiezo = rValue;
        break;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

bool NuTo::LinearPiezoelectric::CheckOutputTypeCompatibility(NuTo::Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
    case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
    case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
    case Constitutive::eOutput::ENGINEERING_STRESS:
    case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
    case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
    case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
    case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
    case Constitutive::eOutput::ELECTRIC_FIELD:
    case Constitutive::eOutput::UPDATE_STATIC_DATA:
    case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}


NuTo::Constitutive::eConstitutiveType NuTo::LinearPiezoelectric::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::LINEAR_PIEZOELECTRIC;
}


void NuTo::LinearPiezoelectric::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Stiffness                     :\n " << this->mStiffness << "\n";
    rLogger << "    Density                       : " << this->mRho << "\n";
    rLogger << "    Piezo                         : " << this->mPiezo << "\n";
    rLogger << "    Permittivity                  : " << this->mPermittivity << "\n";
}


void NuTo::LinearPiezoelectric::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, mRho);
    // Check stiffness tensor components
    // Symmetry test
    for (int ii = 0; ii < 6; ii++)
    {
        for (int jj = 0; jj < ii; jj++)
        {
            if (mStiffness(ii, jj) != mStiffness(jj, ii))
            {
                throw NuTo::Exception(
                        BOOST_CURRENT_FUNCTION,
                        "Stiffness must be symmetric (entry [" + std::to_string(ii) + "," + std::to_string(jj) +
                                "] = " + std::to_string(mStiffness(ii, jj)) + "\n" + "(entry [" + std::to_string(jj) +
                                "," + std::to_string(ii) + "] = " + std::to_string(mStiffness(jj, ii)) + ".");
            }
        }
    }
}
