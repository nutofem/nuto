
#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif // ENABLE_SERIALIZATION

#include "nuto/mechanics/constitutive/laws/PhaseField.h"
#include "nuto/mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataElasticEnergyDensity.h"

#include "nuto/base/Logger.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/constitutive/ConstitutiveBase.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/mechanics/sections/SectionBase.h"
#include "nuto/mechanics/sections/SectionEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/elements/IpDataStaticDataBase.h"

#include "eigen3/Eigen/Eigenvalues"
#include "eigen3/Eigen/Dense"


NuTo::PhaseField::PhaseField(   const double rYoungsModulus,
                                const double rPoissonsRatio,
                                const double rLengthScaleParameter,
                                const double rFractureEnergy,
                                const double rArtificialViscosity)
                                :
                                ConstitutiveBase(),
                                mYoungsModulus                  (rYoungsModulus),
                                mPoissonsRatio                  (rPoissonsRatio),
                                mLengthScaleParameter           (rLengthScaleParameter),
                                mFractureEnergy                 (rFractureEnergy),
                                mArtificialViscosity            (rArtificialViscosity),
                                mLameLambda                     ((rYoungsModulus*rPoissonsRatio) / ((1+rPoissonsRatio)*(1-2*rPoissonsRatio))),
                                mLameMu                         (rYoungsModulus/(2*(1+rPoissonsRatio)))
                                {}


NuTo::ConstitutiveInputMap NuTo::PhaseField::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput, const InterpolationType& rInterpolationType) const
{
    ConstitutiveInputMap constitutiveInputMap;

    constitutiveInputMap[Constitutive::Input::ENGINEERING_STRAIN];
    constitutiveInputMap[Constitutive::Input::CRACK_PHASE_FIELD];

    return constitutiveInputMap;
}

NuTo::Error::eError NuTo::PhaseField::Evaluate1D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}

NuTo::Error::eError NuTo::PhaseField::Evaluate2D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{


    const auto& engineeringStrain   = rConstitutiveInput.at(Constitutive::Input::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& damage              = *rConstitutiveInput.at(Constitutive::Input::CRACK_PHASE_FIELD);

    // calculate coefficients
    double C11, C12, C33;
    switch (rElement->GetSection()->GetType())
    {
    case Section::PLANE_STRAIN:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mYoungsModulus, mPoissonsRatio);
        break;
    case Section::PLANE_STRESS:
        std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(mYoungsModulus, mPoissonsRatio);
        break;
    default:
        throw MechanicsException(std::string("[") + __PRETTY_FUNCTION__ + "[ Invalid type of 2D section behavior found!!!");
    }


    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::Input::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw MechanicsException(__PRETTY_FUNCTION__, "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData = dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);
    int index = calculateStaticData.GetIndexOfPreviousStaticData();
    const ConstitutiveStaticDataElasticEnergyDensity& oldStaticData = *(rElement->GetStaticDataBase(rIp).GetStaticData(index)->AsElasticEnergyDensity());

    // calculate the effective stress
    Eigen::Vector3d effectiveStress;
    effectiveStress[0] = (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
    effectiveStress[1] = (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
    effectiveStress[2] = C33 * engineeringStrain[2];

    // calculate the elastic energy density
    const double elasticEnergyDensity = 0.5*effectiveStress.dot(engineeringStrain);

    ConstitutiveStaticDataElasticEnergyDensity currentStaticData;
    currentStaticData.SetMaxElasticEnergyDensity(std::max(elasticEnergyDensity, oldStaticData.GetMaxElasticEnergyDensity()));

//    Eigen::Matrix2d strainMatrix(2,2);
//    strainMatrix(0,0) = engineeringStrain[0];
//    strainMatrix(0,1) = engineeringStrain[2];

//    strainMatrix(1,0) = engineeringStrain[2];
//    strainMatrix(1,1) = engineeringStrain[1];


//    const double traceStrain = strainMatrix.diagonal().sum();

//    Eigen::Matrix2d strainMatrixPositive(2,2);
//    strainMatrixPositive(0,0) = std::max(engineeringStrain[0],0.0);
//    strainMatrixPositive(0,1) = std::max(engineeringStrain[2],0.0);;
//    strainMatrixPositive(1,0) = std::max(engineeringStrain[2],0.0);;
//    strainMatrixPositive(1,1) = std::max(engineeringStrain[1],0.0);;

//    const double traceStrainPositive = strainMatrixPositive.diagonal().sum();

//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(strainMatrix);

//    //const auto& eigenValues     = eigenSolver.eigenvalues();
//    const auto& eigenVectors    = eigenSolver.eigenvectors();


    bool performUpdateAtEnd = false;
    constexpr double residualEnergyDensity = 1.e-8;

    for (auto itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::Output::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor  = std::pow(1. - damage[0],2) + residualEnergyDensity;

            engineeringStress[0] = factor * effectiveStress[0];
            engineeringStress[1] = factor * effectiveStress[1];
            engineeringStress[2] = factor * effectiveStress[2];

            break;
        }
        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3,3>(itOutput.first, __PRETTY_FUNCTION__);

            // right coefficients are calculated above
            const double factor  = std::pow(1. - damage[0],2) + residualEnergyDensity;

            tangent(0, 0)        = factor * C11;
            tangent(1, 0)        = factor * C12;
            tangent(2, 0)        = 0.;

            tangent(0, 1)        = factor * C12;
            tangent(1, 1)        = factor * C11;
            tangent(2, 1)        = 0.;

            tangent(0, 2)        = 0.;
            tangent(1, 2)        = 0.;
            tangent(2, 2)        = factor * C33;
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor  = std::pow(1. - damage[0],2) + residualEnergyDensity;
            switch (rElement->GetSection()->GetType())
            {
            case Section::PLANE_STRAIN:
                engineeringStress3D[0] = factor * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
                engineeringStress3D[1] = factor * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
                engineeringStress3D[2] = factor *  C12 *(engineeringStrain[0] + engineeringStrain[1]);
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = factor *  C33 * engineeringStrain[2];
                break;
            case Section::PLANE_STRESS:
                engineeringStress3D[0] = factor * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
                engineeringStress3D[1] = factor * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
                engineeringStress3D[2] = 0.;
                engineeringStress3D[3] = 0.;
                engineeringStress3D[4] = 0.;
                engineeringStress3D[5] = factor * C33 * engineeringStrain[2];
                break;

            default:
                throw MechanicsException(__PRETTY_FUNCTION__,"Invalid type of 2D section behavior found!!!");
            }
            break;
        }

        case NuTo::Constitutive::Output::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() = engineeringStrain.As3D(mPoissonsRatio, rElement->GetSection()->GetType());
            break;
        }

        case NuTo::Constitutive::Output::ELASTIC_ENERGY_DAMAGED_PART:
        {
            ConstitutiveIOBase& elasticEnergyLoadTerm = *itOutput.second;

            switch (mPhaseFieldDegradationFunctionType)
            {
            case Constitutive::ePhaseFieldDegradationFunctionType::ISOTROPIC:
            {
                elasticEnergyLoadTerm[0] = currentStaticData.GetMaxElasticEnergyDensity();
                break;
            }
            case Constitutive::ePhaseFieldDegradationFunctionType::ANISOTROPIC_SPECTRAL_DECOMPOSITION:
            {
                //todo: use static data (history variable) for loading/unloading
                //elasticEnergyLoadTerm[0] = 0.5*mLameLambda*std::pow(std::max(traceStrain,0.0),2) + mNu*std::pow(traceStrainPositive,2);
                break;
            }
            default:
                 throw MechanicsException(__PRETTY_FUNCTION__, "Degradation function type not implemented.");
            }

            break;
        }

        case NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_PHASE_FIELD:
        {
            ConstitutiveIOBase& dStressDPhaseField = *itOutput.second;
            dStressDPhaseField.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            const double degradationFactor = 2*(1-damage[0]);

            switch (mPhaseFieldDegradationFunctionType)
            {
            case Constitutive::ePhaseFieldDegradationFunctionType::ISOTROPIC:
            {
                dStressDPhaseField[0] = degradationFactor*effectiveStress[0];
                dStressDPhaseField[1] = degradationFactor*effectiveStress[1];
                dStressDPhaseField[2] = degradationFactor*effectiveStress[2];

                break;
            }
            case Constitutive::ePhaseFieldDegradationFunctionType::ANISOTROPIC_SPECTRAL_DECOMPOSITION:
            {
//                Eigen::Matrix<double,dimension,dimension> engineeringStressPositive;
//                engineeringStressPositive.Zero();

//                for (int i = 0; i < dimension; ++i)
//                {
//                    engineeringStressPositive += eigenVectors.col(i) * eigenVectors.col(i).transpose();
//                }
//                engineeringStressPositive *=mLameLambda*std::max(traceStrain,0.0) + 2*mLameMu*traceStrainPositive;

//                dStressDPhaseField[0] = degradationFactor*engineeringStressPositive.at(0,0);
//                dStressDPhaseField[1] = degradationFactor*engineeringStressPositive.at(1,1);
//                dStressDPhaseField[2] = degradationFactor*engineeringStressPositive.at(1,0);

                break;
            }
            default:
                throw MechanicsException(__PRETTY_FUNCTION__, "Degradation function type not implemented.");
            }

            break;
        }

        case NuTo::Constitutive::Output::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
        {
            // the tangent  is equal to the damaged part of the engineering stress for loading and zero for unloading
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            if (elasticEnergyDensity == currentStaticData.GetMaxElasticEnergyDensity())
            {
                // loading
                switch (mPhaseFieldDegradationFunctionType)
                {
                case Constitutive::ePhaseFieldDegradationFunctionType::ISOTROPIC:
                {
                    tangent[0] = effectiveStress[0];
                    tangent[1] = effectiveStress[1];
                    tangent[2] = effectiveStress[2];

                    break;
                }
                case Constitutive::ePhaseFieldDegradationFunctionType::ANISOTROPIC_SPECTRAL_DECOMPOSITION:
                {

//                    Eigen::Matrix<double,dimension,dimension> engineeringStressPositive;
//                    engineeringStressPositive.Zero();


//                    for (int i = 0; i < dimension; ++i)
//                    {
//                        engineeringStressPositive += eigenVectors.col(i) * eigenVectors.col(i).transpose();
//                    }

//                    engineeringStressPositive *=mLameLambda*std::max(traceStrain,0.0) + 2*mLameMu*traceStrainPositive;

//                    tangent[0] = engineeringStressPositive.at(0,0);
//                    tangent[1] = engineeringStressPositive.at(1,1);
//                    tangent[2] = engineeringStressPositive.at(1,0);
                    break;
                }
                default:
                    throw MechanicsException(__PRETTY_FUNCTION__, "Degradation function type not implemented.");
                }


            } else
            {
                // unloading
                tangent.SetZero();
            }

            break;
        }

        case NuTo::Constitutive::Output::UPDATE_TMP_STATIC_DATA:
        {
            throw MechanicsException(__PRETTY_FUNCTION__,"tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;

        case NuTo::Constitutive::Output::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;

        }
        continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    //update history variables
    if (performUpdateAtEnd)
        rElement->GetStaticData(rIp)->AsElasticEnergyDensity()->SetMaxElasticEnergyDensity(currentStaticData.GetMaxElasticEnergyDensity());

    return Error::SUCCESSFUL;
}

NuTo::Error::eError NuTo::PhaseField::Evaluate3D(ElementBase* rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}


double NuTo::PhaseField::CalculateStaticDataExtrapolationError(ElementBase& rElement, int rIp, const ConstitutiveInputMap& rConstitutiveInput) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}

NuTo::ConstitutiveStaticDataBase* NuTo::PhaseField::AllocateStaticData1D(const ElementBase* rElement) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}

NuTo::ConstitutiveStaticDataBase* NuTo::PhaseField::AllocateStaticData2D(const ElementBase* rElement) const
{
    return new ConstitutiveStaticDataElasticEnergyDensity;
}

NuTo::ConstitutiveStaticDataBase* NuTo::PhaseField::AllocateStaticData3D(const ElementBase* rElement) const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}

bool NuTo::PhaseField::CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::DISPLACEMENTS, Node::DISPLACEMENTS):
        case Node::CombineDofs(Node::DISPLACEMENTS, Node::CRACKPHASEFIELD):
        case Node::CombineDofs(Node::CRACKPHASEFIELD, Node::DISPLACEMENTS):
        case Node::CombineDofs(Node::CRACKPHASEFIELD, Node::CRACKPHASEFIELD):
            return true;
        default:
            return false;
        }
        break;
    case 1:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::CRACKPHASEFIELD, Node::CRACKPHASEFIELD):
            return true;
        default:
            return false;
        }
        break;
    default:
        return false;
    }

}

double NuTo::PhaseField::GetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier) const
{
    switch(rIdentifier)
    {
    case Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY:
        return mArtificialViscosity;
    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
        return mFractureEnergy;
    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
        return mPoissonsRatio;
    case Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER:
        return mLengthScaleParameter;
    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
        return mYoungsModulus;

    default:
        throw MechanicsException(__PRETTY_FUNCTION__,"Constitutive law does not have the requested variable");
    }
}

void NuTo::PhaseField::SetParameterDouble(Constitutive::eConstitutiveParameter rIdentifier, double rValue)
{
//    switch (rIdentifier)
//    {
//    case Constitutive::eConstitutiveParameter::ARTIFICIAL_VISCOSITY:
//        mArtificialViscosity = rValue;
//        break;
//    case Constitutive::eConstitutiveParameter::FRACTURE_ENERGY:
//        mFractureEnergy = rValue;
//        break;
//    case Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER:
//        mLengthScaleParameter = rValue;
//        break;
//    case Constitutive::eConstitutiveParameter::POISSONS_RATIO:
//        mPoissonsRatio = rValue;
//        break;
//    case Constitutive::eConstitutiveParameter::YOUNGS_MODULUS:
//        mYoungsModulus = rValue;
//        break;
//    default:
//        throw MechanicsException(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
//    }
//    SetParametersValid();
}

NuTo::Constitutive::eConstitutiveType NuTo::PhaseField::GetType() const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented.");
}

void NuTo::PhaseField::CheckParameters() const
{
}

bool NuTo::PhaseField::CheckElementCompatibility(Element::eElementType rElementType) const
{
    switch (rElementType)
    {
    case NuTo::Element::CONTINUUMELEMENT:
        return true;
    default:
        return false;
    }
}

void NuTo::PhaseField::Info(unsigned short rVerboseLevel, Logger& rLogger) const
{
  this->ConstitutiveBase::Info(rVerboseLevel, rLogger);
  rLogger << "    Young's modulus         : " << mYoungsModulus << "\n";
  rLogger << "    Poisson's ratio         : " << mPoissonsRatio << "\n";
  rLogger << "    fracture energy         : " << mFractureEnergy << "\n";
  rLogger << "    artificial viscosity    : " << mArtificialViscosity << "\n";
  rLogger << "    length scale parameter  : " << mLengthScaleParameter << "\n";
}
