#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/staticData/IPConstitutiveLawWithoutData.h"

#include "mechanics/constitutive/laws/LinearPiezoelectric.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "base/Logger.h"
#include "mechanics/MechanicsException.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutivePlaneState.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

#include <boost/current_function.hpp>

NuTo::LinearPiezoelectric::LinearPiezoelectric() :
        ConstitutiveBase()
{
    mStiffness << 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0;
    mRho = 0.;

    mPermittivity << 0, 0, 0,
                     0, 0, 0,
                     0, 0, 0;

    mPiezo <<  0, 0, 0,
               0, 0, 0,
               0, 0, 0,
               0, 0, 0,
               0, 0, 0,
               0, 0, 0;

    SetParametersValid();
}

//#ifdef ENABLE_SERIALIZATION
////! @brief serializes the class
////! @param ar         archive
////! @param version    version
//template<class Archive>
//void NuTo::LinearPiezoelectric::serialize(Archive & ar, const unsigned int version)
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
//BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LinearPiezoelectric)
//#endif // ENABLE_SERIALIZATION

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::LinearPiezoelectric::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<LinearPiezoelectric>>(*this);
}


NuTo::ConstitutiveInputMap NuTo::LinearPiezoelectric::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
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
//            throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] output object " + Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

namespace NuTo // template specialization in same namespace as definition
{
template<>
void NuTo::LinearPiezoelectric::Evaluate<1>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    std::cout << "Not implemented" << std::endl;
}


template<>
void NuTo::LinearPiezoelectric::Evaluate<2>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    std::cout << "Not implemented" << std::endl;
}


template<>
void NuTo::LinearPiezoelectric::Evaluate<3>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    InputData<3> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch (itInput.first)
        {
        case Constitutive::eInput::ELECTRIC_FIELD:
            inputData.mElectricField = static_cast<ConstitutiveVector<3>*>(itInput.second.get())->AsVector();
            break;
        default:
            continue;
        }
    }


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

                Eigen::Matrix<double, 3, 1>& electricField = *static_cast<ConstitutiveVector<3>*>(itOutput.second.get());
                electricField = inputData.mElectricField;

                // transform vectors to eigen vectors
                // let eigen do the matrix multiply
                // and transform back: This can certainly be done better.
                Eigen::VectorXd engStress(6);
                Eigen::VectorXd engStrain(6);
                for (int ii=0; ii<6; ii++) {
                    engStrain(ii) = engineeringStrain[ii];
                }

                engStress = this->mStiffness * engStrain - (this->mPiezo).transpose() * electricField;

                for (int ii=0; ii<6; ii++) {
                    engineeringStress[ii] = engStress(ii);
                }

                break;
            }
        case NuTo::Constitutive::eOutput::ELECTRIC_FIELD:
        {
            Eigen::Matrix<double, 3, 1>& electricField = *static_cast<ConstitutiveVector<3>*>(itOutput.second.get());
            electricField = inputData.mElectricField;
            break;
        }
        case NuTo::Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            ConstitutiveIOBase& electricDisplacement = *itOutput.second;
            electricDisplacement.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain =
                rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();

            Eigen::Matrix<double, 3, 1>& electricField = *static_cast<ConstitutiveVector<3>*>(itOutput.second.get());
            electricField = inputData.mElectricField;

            // transform vectors to eigen vectors
            // let eigen do the matrix multiply
            // and transform back: This can certainly be done better.
            Eigen::VectorXd engStrain(6);
            Eigen::VectorXd elDispl(3);
            for (int ii=0; ii<6; ii++) {
                engStrain(ii) = engineeringStrain[ii];
            }

            elDispl = this->mPiezo * engStrain + this->mPermittivity * electricField;

            for (int ii=0; ii<3; ii++) {
                electricDisplacement[ii] = elDispl(ii);
            }

            break;
        }
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<6, 6>(itOutput.first, __PRETTY_FUNCTION__);

                for (int ii=0; ii<6; ii++) {
                    for (int jj=0; jj<6; jj++) {
                        tangent(ii,jj) = this->mStiffness(ii,jj);
                    }
                }

                break;
            }
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ELECTRIC_FIELD:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<6, 3>(itOutput.first, __PRETTY_FUNCTION__);

                for (int ii=0; ii<6; ii++) {
                    for (int jj=0; jj<3; jj++) {
                        tangent(ii,jj) = ((this->mPiezo).transpose())(ii,jj);
                    }
                }

                break;
            }
            case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ENGINEERING_STRAIN:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<3, 6>(itOutput.first, __PRETTY_FUNCTION__);

                for (int ii=0; ii<3; ii++) {
                    for (int jj=0; jj<6; jj++) {
                        tangent(ii,jj) = this->mPiezo(ii,jj);
                    }
                }

                break;
            }
            case NuTo::Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

                for (int ii=0; ii<3; ii++) {
                    for (int jj=0; jj<3; jj++) {
                        tangent(ii,jj) = this->mPermittivity(ii,jj);
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
                //nothing to be done for update routine
                continue;
            }
            default:
                continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

} // namespace NuTo


bool NuTo::LinearPiezoelectric::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow,
                                                                          NuTo::Node::eDof rDofCol,
                                                                          int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
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
    switch(rIdentifier)
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
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}


void NuTo::LinearPiezoelectric::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        this->mRho = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

Eigen::VectorXd NuTo::LinearPiezoelectric::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::STIFFNESS:
    {
        Eigen::VectorXd stiffnessFlattened(36);
        int count = 0;
        for (int ii=0; ii< 6; ii++) {
            for (int jj=0; jj< 6; jj++) {
                stiffnessFlattened(count) = this->mStiffness(ii,jj);
                count++;
            }
        }
        return stiffnessFlattened;
    }
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        Eigen::VectorXd dielectricTensorFlattened(9);
        int count = 0;
        for (int ii=0; ii< 3; ii++) {
            for (int jj=0; jj< 3; jj++) {
                dielectricTensorFlattened(count) = this->mPermittivity(ii,jj);
                count++;
            }
        }
        return dielectricTensorFlattened;
    }
    case Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR:
    {
        Eigen::VectorXd piezoelectricTensorFlattened(18);
        int count = 0;
        for (int ii=0; ii< 3; ii++) {
            for (int jj=0; jj< 6; jj++) {
                piezoelectricTensorFlattened(count) = this->mPiezo(ii,jj);
                count++;
            }
        }
        return piezoelectricTensorFlattened;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::LinearPiezoelectric::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue)
{
    //ConstitutiveBase::CheckParameterFullVector(rIdentifier, rValue);
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::STIFFNESS:
    {
        int count = 0;
        for (int ii=0; ii< 6; ii++) {
            for (int jj=0; jj< 6; jj++) {
                this->mStiffness(ii,jj) = rValue(count);
                count++;
            }
        }
        break;
    }
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        int count = 0;
        for (int ii=0; ii< 3; ii++) {
            for (int jj=0; jj< 3; jj++) {
                this->mPermittivity(ii,jj) = rValue(count);
                count++;
            }
        }
        break;
    }
    case Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR:
    {
        int count = 0;
        for (int ii=0; ii< 3; ii++) {
            for (int jj=0; jj< 6; jj++) {
                this->mPiezo(ii,jj) = rValue(count);
                count++;
            }
        }
        break;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
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
//    case Constitutive::eOutput::UPDATE_STATIC_DATA:
//    case Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
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
    for (int ii=0; ii<6; ii++) {
        for (int jj=0; jj<ii; jj++) {
            if (mStiffness(ii,jj) != mStiffness(jj,ii)) {
                throw NuTo::MechanicsException(BOOST_CURRENT_FUNCTION , "Stiffness must be symmetric (entry ["
                                               + std::to_string(ii) + "," + std::to_string(jj) + "] = "
                                               + std::to_string(mStiffness(ii,jj)) + "\n"
                                               + "(entry [" + std::to_string(jj) + "," + std::to_string(ii) + "] = "
                                               + std::to_string(mStiffness(jj,ii)) + ".");
            }
        }
    }
}
