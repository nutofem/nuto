// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"

#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"


NuTo::LinearElasticEngineeringStress::LinearElasticEngineeringStress() :
        ConstitutiveBase()
{
    mE = 0.;
    mNu = 0.;
    mRho = 0.;
    mThermalExpansionCoefficient = 0.;
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



NuTo::ConstitutiveInputMap NuTo::LinearElasticEngineeringStress::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            constitutiveInputMap[Constitutive::Input::ENGINEERING_STRAIN];

            if (rInterpolationType.IsConstitutiveInput(Node::TEMPERATURE))
                constitutiveInputMap[Constitutive::Input::TEMPERATURE];

            break;
        }
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
            constitutiveInputMap[Constitutive::Input::ENGINEERING_STRAIN];
            break;
        // no inputs needed for:
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
            break;
        default:
            continue;
//            ProcessUnhandledOutput(__PRETTY_FUNCTION__,itOutput.first);
//            throw MechanicsException(std::string("[")+__PRETTY_FUNCTION__+"] output object " + Constitutive::OutputToString(itOutput.first) + " cannot be calculated by this constitutive law.");
        }
    }

    return constitutiveInputMap;
}

NuTo::Error::eError NuTo::LinearElasticEngineeringStress::Evaluate1D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain1D();
            auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<1>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<1>(itOutput.first, __PRETTY_FUNCTION__);

            engineeringStress[0] = mE * elasticEngineeringStrain[0];
            break;
        }
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain1D();
            auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<1>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

            engineeringStress3D.SetZero();
            engineeringStress3D[0] = mE * elasticEngineeringStrain[0];
            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<1,1>(itOutput.first, __PRETTY_FUNCTION__);

            tangent(0,0) = mE;
            break;
        }
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain1D().As3D(mNu);
            break;
        }
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            //nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    //update history variables but for linear elastic, there is nothing to do
    return Error::SUCCESSFUL;
}





NuTo::Error::eError NuTo::LinearElasticEngineeringStress::Evaluate2D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    assert(rElement->GetSection() != nullptr); // section needed to determine plane_stress/plane_strain

    // calculate coefficients
    double C11, C12, C33;
    switch (rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
        break;
    case Section::PLANE_STRESS:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(mE, mNu);
        break;
    default:
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "[ Invalid type of 2D section behavior found!!!");
    }


    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
            auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<2>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

            engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
            engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
            engineeringStress[2] = C33 * elasticEngineeringStrain[2];
            break;
        }
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE: // for visualization
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
            auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<2>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
            {
                engineeringStress3D[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
                engineeringStress3D[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
                engineeringStress3D[2] = C12 * (elasticEngineeringStrain[0] + elasticEngineeringStrain[1]);
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = C33 * elasticEngineeringStrain[2];
                break;
            }
            case Section::PLANE_STRESS:
            {
                engineeringStress3D[0] = C11 * elasticEngineeringStrain[0] + C12 * elasticEngineeringStrain[1];
                engineeringStress3D[1] = C11 * elasticEngineeringStrain[1] + C12 * elasticEngineeringStrain[0];
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = C33 * elasticEngineeringStrain[2];
                break;
            }
            default:
                throw MechanicsException(std::string("[") +__PRETTY_FUNCTION__ + "[ Invalid type of 2D section behavior found!!!");
            }

            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
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
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE: // for visualization
        {
            itOutput.second->AsEngineeringStrain3D() = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain2D().As3D(mNu, rElement->GetSection()->GetType());
        }
            break;
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            //nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
    //update history variables but for linear elastic, there is nothing to do
    return Error::SUCCESSFUL;
}





NuTo::Error::eError NuTo::LinearElasticEngineeringStress::Evaluate3D(
        ElementBase* rElement, int rIp,
        const ConstitutiveInputMap& rConstitutiveInput,
        const ConstitutiveOutputMap& rConstitutiveOutput)
{
    double C11=0.0, C12=0.0, C44=0.0;
    if (rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS) != rConstitutiveOutput.end()
     || rConstitutiveOutput.find(NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN) != rConstitutiveOutput.end()
     || rConstitutiveOutput.find(NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE) != rConstitutiveOutput.end())
    {
        std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(mE, mNu);
    }

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const auto& engineeringStrain = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain3D();
            auto elasticEngineeringStrain = EngineeringStressHelper::CalculateElasticEngineeringStrain<3>(engineeringStrain, *rElement->GetInterpolationType(), rConstitutiveInput, mThermalExpansionCoefficient);

            engineeringStress[0] = C11 * elasticEngineeringStrain[0] + C12 * (elasticEngineeringStrain[1] + elasticEngineeringStrain[2]);
            engineeringStress[1] = C11 * elasticEngineeringStrain[1] + C12 * (elasticEngineeringStrain[0] + elasticEngineeringStrain[2]);
            engineeringStress[2] = C11 * elasticEngineeringStrain[2] + C12 * (elasticEngineeringStrain[0] + elasticEngineeringStrain[1]);
            engineeringStress[3] = C44 * elasticEngineeringStrain[3];
            engineeringStress[4] = C44 * elasticEngineeringStrain[4];
            engineeringStress[5] = C44 * elasticEngineeringStrain[5];
            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<6,6>(itOutput.first, __PRETTY_FUNCTION__);

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
        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain3D();
            break;
        }
        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            //nothing to be done for update routine
            continue;
        }
        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }

    //update history variables but for linear elastic, there is nothing to do

    return Error::SUCCESSFUL;
}




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
        case Node::CombineDofs(Node::DISPLACEMENTS, Node::DISPLACEMENTS):
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
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        return this->mThermalExpansionCoefficient;
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
    case Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT:
        this->mThermalExpansionCoefficient = rValue;
        break;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        this->mE = rValue;
        break;
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}





bool NuTo::LinearElasticEngineeringStress::CheckOutputTypeCompatibility(NuTo::Constitutive::Output::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
    case Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
    case Constitutive::Output::ENGINEERING_STRESS:
    case Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
    case Constitutive::Output::UPDATE_STATIC_DATA:
    case Constitutive::Output::UPDATE_TMP_STATIC_DATA:
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
    return NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS;
}



bool NuTo::LinearElasticEngineeringStress::CheckElementCompatibility(NuTo::Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::CONTINUUMELEMENT:
    case NuTo::Element::CONTINUUMELEMENTIGA:
    case NuTo::Element::CONTINUUMBOUNDARYELEMENTCONSTRAINEDCONTROLNODE:
    case NuTo::Element::ELEMENT1DINXD:
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
    rLogger << "    thermal expansion coefficient : " << this->mThermalExpansionCoefficient << "\n";
}



void NuTo::LinearElasticEngineeringStress::CheckParameters() const
{
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, mE);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, mNu);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, mRho);
    ConstitutiveBase::CheckParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, mThermalExpansionCoefficient);
}



