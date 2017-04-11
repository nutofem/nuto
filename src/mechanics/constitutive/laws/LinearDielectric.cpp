#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/laws/LinearDielectric.h"
#include "base/Logger.h"
#include "mechanics/MechanicsException.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveVector.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include "mechanics/elements/ElementBase.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/nodes/NodeEnum.h"

using namespace NuTo;

LinearDielectric::LinearDielectric() : ConstitutiveBase()
{
//    mPermittivity << 1, 0, 0,
//                     0, 1, 0,
//                     0, 0, 1;
    SetParametersValid();
}

#ifdef ENABLE_SERIALIZATION
template <class Archive> void LinearDielectric::serialize(Archive& ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize LinearDielectric" << std::endl;
#endif
    ar& BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstitutiveBase) & BOOST_SERIALIZATION_NVP(mK) &
            BOOST_SERIALIZATION_NVP(mCt) & BOOST_SERIALIZATION_NVP(mRho);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize LinearDielectric" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(LinearDielectric)
#endif // ENABLE_SERIALIZATION

ConstitutiveInputMap LinearDielectric::GetConstitutiveInputs(
        const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        case Constitutive::eOutput::ELECTRIC_FIELD:
        {
            constitutiveInputMap[Constitutive::eInput::ELECTRIC_FIELD];
            break;
        }
        case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            // no input needed for the linear case
            // the result is just a constant
            break;
        }
        default:
            continue;
        }
    }

    return constitutiveInputMap;
}

bool LinearDielectric::CheckDofCombinationComputable(Node::eDof dofRow, Node::eDof dofCol, int timeDerivative) const
{
    if ((timeDerivative == 0) && (dofRow == Node::eDof::ELECTRICPOTENTIAL and dofCol == Node::eDof::ELECTRICPOTENTIAL))
    {
            return true;
    }
    else return false;
}

template <int TDim>
void LinearDielectric::Evaluate(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    Eigen::Matrix<double,3,TDim> cropMatrix;
    for(int ii=0; ii<3; ii++) {
        for (int jj=0; jj <TDim; jj++) {
            if (ii == jj) {
                cropMatrix(ii,jj) = 1;
            } else {
                cropMatrix(ii,jj) = 0;
            }
        }
    }
    auto perm = cropMatrix.transpose() * mPermittivity * cropMatrix;

    InputData<TDim> inputData;
    for (auto& itInput : rConstitutiveInput)
    {
        switch (itInput.first)
        {
        case Constitutive::eInput::ELECTRIC_FIELD:
            inputData.mElectricField = static_cast<ConstitutiveVector<TDim>*>(itInput.second.get())->AsVector();
            break;
        default:
            continue;
        }
    }

    for (const auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case Constitutive::eOutput::ELECTRIC_FIELD:
        {
            Eigen::Matrix<double, TDim, 1>& electricField = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            electricField = inputData.mElectricField;
            break;
        }

        case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
        {
            Eigen::Matrix<double, TDim, 1>& electricDisplacement = *static_cast<ConstitutiveVector<TDim>*>(itOutput.second.get());
            electricDisplacement = perm * inputData.mElectricField;
            break;
        }
        case Constitutive::eOutput::D_ELECTRIC_DISPLACEMENT_D_ELECTRIC_FIELD:
        {
            Eigen::Matrix<double, TDim, TDim>& permittivity =
                    *static_cast<ConstitutiveMatrix<TDim, TDim>*>(itOutput.second.get());
            permittivity = perm;
            break;
        }

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }
}

void NuTo::LinearDielectric::CheckParameters() const
{

}

bool LinearDielectric::CheckHaveParameter(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch (rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

Eigen::VectorXd NuTo::LinearDielectric::GetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR:
    {
        Eigen::VectorXd permittivityFlattened(9);
        int count = 0;
        for (int ii=0; ii< 3; ii++) {
            for (int jj=0; jj< 3; jj++) {
                permittivityFlattened(count) = this->mPermittivity(ii,jj);
                count++;
            }
        }
        return permittivityFlattened;
    }
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::LinearDielectric::SetParameterFullVectorDouble(NuTo::Constitutive::eConstitutiveParameter rIdentifier, Eigen::VectorXd rValue)
{
    //ConstitutiveBase::CheckParameterFullVector(rIdentifier, rValue);
    switch(rIdentifier)
    {
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
    default:
        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

bool LinearDielectric::CheckOutputTypeCompatibility(Constitutive::eOutput rOutputEnum) const
{
    switch (rOutputEnum)
    {
    case Constitutive::eOutput::ELECTRIC_DISPLACEMENT:
    {
        return true;
    }
    default:
    {
        return false;
    }
    }
}

Constitutive::eConstitutiveType LinearDielectric::GetType() const
{
    return Constitutive::eConstitutiveType::LINEAR_DIELECTRIC;
}

void LinearDielectric::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
    this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
    rLogger << "    Relative electric permittivity          : " << this->mPermittivity << "\n";
}

template void LinearDielectric::Evaluate<1>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearDielectric::Evaluate<2>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
template void LinearDielectric::Evaluate<3>(
        const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput);
