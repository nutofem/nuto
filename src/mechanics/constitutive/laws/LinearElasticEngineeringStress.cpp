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

#include "mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
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
#include "mechanics/sections/SectionBase.h"
#include "mechanics/sections/SectionEnum.h"


NuTo::LinearElasticEngineeringStress::LinearElasticEngineeringStress() :
        ConstitutiveBase()
{
    mE = 0.;
    mNu = 0.;
    mRho = 0.;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
//! @brief serializes the class
//! @param ar         archive
//! @param version    version
template<class Archive>
void NuTo::LinearElasticEngineeringStress::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LinearElasticEngineeringStress" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase)
    & BOOST_SERIALIZATION_NVP(mE)
    & BOOST_SERIALIZATION_NVP(mNu)
    & BOOST_SERIALIZATION_NVP(mRho)
    & BOOST_SERIALIZATION_NVP(mThermalExpansionCoefficient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LinearElasticEngineeringStress" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::LinearElasticEngineeringStress)
#endif // ENABLE_SERIALIZATION

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::LinearElasticEngineeringStress::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLawWithoutData<LinearElasticEngineeringStress>>(*this);
}


NuTo::ConstitutiveInputMap NuTo::LinearElasticEngineeringStress::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
            break;
        // no inputs needed for:
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
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
void NuTo::LinearElasticEngineeringStress::Evaluate<1>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
            case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
            {
                const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

                ConstitutiveIOBase& engineeringStress = *itOutput.second;
                engineeringStress.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

                engineeringStress[0] = mE * engineeringStrain[0];
                break;
            }
            case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
            {
                ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
                engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

                const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D();

                engineeringStress3D.SetZero();
                engineeringStress3D[0] = mE * engineeringStrain[0];
                break;
            }
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<1, 1>(itOutput.first, __PRETTY_FUNCTION__);

                tangent(0, 0) = mE;
                break;
            }
            case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
            {
                itOutput.second->AsEngineeringStrain3D() =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain1D().As3D(mNu);
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


template<>
void NuTo::LinearElasticEngineeringStress::Evaluate<2>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    const auto& planeState =
        *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());

    // calculate coefficients
    double C11, C12, C33;
    switch (planeState.GetPlaneState())
    {
        case ePlaneState::PLANE_STRAIN:
            std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
            break;
        case ePlaneState::PLANE_STRESS:
            std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(mE, mNu);
            break;
        default:
            throw MechanicsException(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");
    }


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

                engineeringStress[0] = C11 * engineeringStrain[0] + C12 * engineeringStrain[1];
                engineeringStress[1] = C11 * engineeringStrain[1] + C12 * engineeringStrain[0];
                engineeringStress[2] = C33 * engineeringStrain[2];
                break;
            }
            case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE: // for visualization
            {
                ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
                engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

                const auto& engineeringStrain =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();

                switch (planeState.GetPlaneState())
                {
                    case ePlaneState::PLANE_STRAIN:
                    {
                        engineeringStress3D[0] = C11 * engineeringStrain[0] + C12 * engineeringStrain[1];
                        engineeringStress3D[1] = C11 * engineeringStrain[1] + C12 * engineeringStrain[0];
                        engineeringStress3D[2] = C12 * (engineeringStrain[0] + engineeringStrain[1]);
                        engineeringStress3D[3] = 0.;
                        engineeringStress3D[4] = 0.;
                        engineeringStress3D[5] = C33 * engineeringStrain[2];
                        break;
                    }
                    case ePlaneState::PLANE_STRESS:
                    {
                        engineeringStress3D[0] = C11 * engineeringStrain[0] + C12 * engineeringStrain[1];
                        engineeringStress3D[1] = C11 * engineeringStrain[1] + C12 * engineeringStrain[0];
                        engineeringStress3D[2] = 0.;
                        engineeringStress3D[3] = 0.;
                        engineeringStress3D[4] = 0.;
                        engineeringStress3D[5] = C33 * engineeringStrain[2];
                        break;
                    }
                    default:
                        throw MechanicsException(
                            std::string("[") + __PRETTY_FUNCTION__ + "[ Invalid type of 2D section behavior found!!!");
                }

                break;
            }
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

                tangent(0, 0) = C11;
                tangent(1, 0) = C12;
                tangent(2, 0) = 0;

                tangent(0, 1) = C12;
                tangent(1, 1) = C11;
                tangent(2, 1) = 0;

                tangent(0, 2) = 0.;
                tangent(1, 2) = 0.;
                tangent(2, 2) = C33;
                break;
            }
            case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE: // for visualization
            {
                itOutput.second->AsEngineeringStrain3D() =
                    rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D()
                        .As3D(mNu, planeState.GetPlaneState());
            }
                break;
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


template<>
void NuTo::LinearElasticEngineeringStress::Evaluate<3>(
    const ConstitutiveInputMap& rConstitutiveInput,
    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    double C11 = 0.0, C12 = 0.0, C44 = 0.0;
    if (rConstitutiveOutput.find(NuTo::Constitutive::eOutput::ENGINEERING_STRESS) != rConstitutiveOutput.end()
        || rConstitutiveOutput.find(NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN)
            != rConstitutiveOutput.end()
        || rConstitutiveOutput.find(NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE)
            != rConstitutiveOutput.end())
    {
        std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
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

                engineeringStress[0] = C11 * engineeringStrain[0] + C12 * (engineeringStrain[1] + engineeringStrain[2]);
                engineeringStress[1] = C11 * engineeringStrain[1] + C12 * (engineeringStrain[0] + engineeringStrain[2]);
                engineeringStress[2] = C11 * engineeringStrain[2] + C12 * (engineeringStrain[0] + engineeringStrain[1]);
                engineeringStress[3] = C44 * engineeringStrain[3];
                engineeringStress[4] = C44 * engineeringStrain[4];
                engineeringStress[5] = C44 * engineeringStrain[5];
                break;
            }
            case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
            {
                ConstitutiveIOBase& tangent = *itOutput.second;
                tangent.AssertIsMatrix<6, 6>(itOutput.first, __PRETTY_FUNCTION__);

                tangent.SetZero();

                // C11 diagonal:
                tangent(0, 0) = C11;
                tangent(1, 1) = C11;
                tangent(2, 2) = C11;

                // C12 off diagonals:
                tangent(0, 1) = C12;
                tangent(0, 2) = C12;
                tangent(1, 0) = C12;
                tangent(1, 2) = C12;
                tangent(2, 0) = C12;
                tangent(2, 1) = C12;

                // C44 diagonal:
                tangent(3, 3) = C44;
                tangent(4, 4) = C44;
                tangent(5, 5) = C44;

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


bool NuTo::LinearElasticEngineeringStress::CheckDofCombinationComputable(NuTo::Node::eDof rDofRow,
                                                                          NuTo::Node::eDof rDofCol,
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


bool NuTo::LinearElasticEngineeringStress::CheckHaveParameter(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
        case Constitutive::eConstitutiveParameter::DENSITY:
        case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        {
            return true;
        }
        default:
        {
            return false;
        }
    }
}


double NuTo::LinearElasticEngineeringStress::GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        return this->mRho;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return this->mNu;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return this->mE;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}


void NuTo::LinearElasticEngineeringStress::SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
    ConstitutiveBase::CheckParameterDouble(rIdentifier, rValue);
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DENSITY:
        this->mRho = rValue;
        break;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        this->mNu = rValue;
        break;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        this->mE = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}


bool NuTo::LinearElasticEngineeringStress::CheckOutputTypeCompatibility(NuTo::Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
    case Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
    case Constitutive::eOutput::ENGINEERING_STRESS:
    case Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
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


NuTo::Constitutive::eConstitutiveType NuTo::LinearElasticEngineeringStress::GetType() const
{
    return NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS;
}


bool NuTo::LinearElasticEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::eElementType::CONTINUUMELEMENT:
    case NuTo::Element::eElementType::CONTINUUMELEMENTIGA:
    case NuTo::Element::eElementType::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
    case NuTo::Element::eElementType::ELEMENT1DINXD:
        return true;
    default:
        return false;
    }
}


void NuTo::LinearElasticEngineeringStress::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Young's modulus               : " << this->mE << "\n";
    rLogger << "    Poisson's ratio               : " << this->mNu << "\n";
    rLogger << "    Density                       : " << this->mRho << "\n";
}


void NuTo::LinearElasticEngineeringStress::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, mE);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, mNu);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, mRho);
}
