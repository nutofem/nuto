#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"

#include "base/Logger.h"
#include "base/Exception.h"
#include "mechanics/constitutive/ConstitutiveBase.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#include "mechanics/nodes/NodeEnum.h"

#include <eigen3/Eigen/Core>

#include "mechanics/constitutive/staticData/IPConstitutiveLaw.h"


NuTo::PhaseField::PhaseField(const double rYoungsModulus, const double rPoissonsRatio,
                             const double rLengthScaleParameter, const double rFractureEnergy,
                             const double rArtificialViscosity,
                             const Constitutive::ePhaseFieldEnergyDecomposition rEnergyDecomposition)
    : ConstitutiveBase()
    , mYoungsModulus(rYoungsModulus)
    , mPoissonsRatio(rPoissonsRatio)
    , mLengthScaleParameter(rLengthScaleParameter)
    , mFractureEnergy(rFractureEnergy)
    , mArtificialViscosity(rArtificialViscosity)
    , mLameLambda((rYoungsModulus * rPoissonsRatio) / ((1 + rPoissonsRatio) * (1 - 2 * rPoissonsRatio)))
    , mLameMu(rYoungsModulus / (2 * (1 + rPoissonsRatio)))
    , mEnergyDecomposition(rEnergyDecomposition)
{
    assert(mYoungsModulus > 0.0);
    assert(mLengthScaleParameter > 0.0);
    assert(mFractureEnergy > 0.0);
    assert(mArtificialViscosity >= 0.0);
    assert(mPoissonsRatio >= 0.0 and mPoissonsRatio < 0.5);
}

std::unique_ptr<NuTo::Constitutive::IPConstitutiveLawBase> NuTo::PhaseField::CreateIPLaw()
{
    return std::make_unique<Constitutive::IPConstitutiveLaw<PhaseField>>(*this, 0.);
}

NuTo::ConstitutiveInputMap NuTo::PhaseField::GetConstitutiveInputs(const ConstitutiveOutputMap&) const
{
    ConstitutiveInputMap constitutiveInputMap;

    constitutiveInputMap[Constitutive::eInput::ENGINEERING_STRAIN];
    constitutiveInputMap[Constitutive::eInput::CRACK_PHASE_FIELD];

    return constitutiveInputMap;
}


// template this function over dimension and reduce the number of input arguments once its working
double NuTo::PhaseField::CalculateComponentsSpectralDecompositionDStressDStrain(
        int rI, int rJ, int rK, int rL, const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 2, 2>>& rEigenSolver,
        std::function<double(double, double)> rRampFunction, std::function<bool(double)> rStepFunction)
{
    const int i = rI;
    const int j = rJ;
    const int k = rK;
    const int l = rL;

    Eigen::Matrix2i kroneckerDelta = Eigen::Matrix2i::Identity();

    const Eigen::Matrix<double, 2, 1> eigenVec00 = rEigenSolver.eigenvectors().col(0);
    const Eigen::Matrix<double, 2, 1> eigenVec01 = rEigenSolver.eigenvectors().col(1);


    const double eigenVal00 = rEigenSolver.eigenvalues()(0, 0);
    const double eigenVal01 = rEigenSolver.eigenvalues()(1, 0);
    const double traceStrain = eigenVal00 + eigenVal01;


    double Cijkl =
            mLameLambda * rStepFunction(traceStrain) * kroneckerDelta(i, j) * kroneckerDelta(k, l) +
            2 * mLameMu * rStepFunction(eigenVal00) * eigenVec00[i] * eigenVec00[j] * eigenVec00[k] * eigenVec00[l] +
            2 * mLameMu * rStepFunction(eigenVal01) * eigenVec01[i] * eigenVec01[j] * eigenVec01[k] * eigenVec01[l];

    if (eigenVal00 != eigenVal01)
    {
        Cijkl += mLameMu * (

                                   rRampFunction(eigenVal00, 0.0) / (eigenVal00 - eigenVal01) *
                                           (

                                                   eigenVec00[j] * eigenVec00[l] * eigenVec01[i] * eigenVec01[k]

                                                   + eigenVec00[j] * eigenVec00[k] * eigenVec01[i] * eigenVec01[l]

                                                   + eigenVec00[i] * eigenVec00[l] * eigenVec01[j] * eigenVec01[k]

                                                   + eigenVec00[i] * eigenVec00[k] * eigenVec01[j] * eigenVec01[l]

                                                   )

                                   +
                                   rRampFunction(eigenVal01, 0.0) / (eigenVal01 - eigenVal00) *
                                           (

                                                   eigenVec01[j] * eigenVec01[l] * eigenVec00[i] * eigenVec00[k]

                                                   + eigenVec01[j] * eigenVec01[k] * eigenVec00[i] * eigenVec00[l]

                                                   + eigenVec01[i] * eigenVec01[l] * eigenVec00[j] * eigenVec00[k]

                                                   + eigenVec01[i] * eigenVec01[k] * eigenVec00[j] * eigenVec00[l]

                                                   )

                                           );
    }

    return Cijkl;
}

void NuTo::PhaseField::CalculateSpectralDecompositionDStressDStrain(
        ConstitutiveIOBase& tangent, const double factor,
        const Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d>& eigenSolver)
{
    auto positiveRampFunction = [](double a, double b) { return std::max<double>(a, b); };
    auto negativeRampFunction = [](double a, double b) { return std::min<double>(a, b); };

    auto positiveStepFunction = [](double a) { return (a >= 0.0); };
    auto negativeStepFunction = [](double a) { return (a < 0.0); };


    int voigtI = 0;
    int voigtJ = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);

    voigtI = 1;
    voigtJ = 1;
    i = 1;
    j = 1;
    k = 1;
    l = 1;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);

    voigtI = 2;
    voigtJ = 2;
    i = 0;
    j = 1;
    k = 0;
    l = 1;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);


    voigtI = 0;
    voigtJ = 1;
    i = 0;
    j = 0;
    k = 1;
    l = 1;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);

    tangent(voigtJ, voigtI) = tangent(voigtI, voigtJ);

    voigtI = 0;
    voigtJ = 2;
    i = 0;
    j = 0;
    k = 0;
    l = 1;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);

    tangent(voigtJ, voigtI) = tangent(voigtI, voigtJ);

    voigtI = 1;
    voigtJ = 2;
    i = 1;
    j = 1;
    k = 0;
    l = 1;
    tangent(voigtI, voigtJ) = factor * CalculateComponentsSpectralDecompositionDStressDStrain(
                                               i, j, k, l, eigenSolver, positiveRampFunction, positiveStepFunction);
    tangent(voigtI, voigtJ) += CalculateComponentsSpectralDecompositionDStressDStrain(
            i, j, k, l, eigenSolver, negativeRampFunction, negativeStepFunction);

    tangent(voigtJ, voigtI) = tangent(voigtI, voigtJ);
}


double NuTo::PhaseField::Evaluate2DIsotropic(const double oldEnergyDensity,
                                             const ConstitutiveInputMap& rConstitutiveInput,
                                             const ConstitutiveOutputMap& rConstitutiveOutput)
{

    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& damage = *rConstitutiveInput.at(Constitutive::eInput::CRACK_PHASE_FIELD);

    // calculate coefficients
    double C11, C12, C33;
    std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients3D(mYoungsModulus, mPoissonsRatio);

    // calculate the effective stress
    Eigen::Vector3d effectiveStress;
    effectiveStress[0] = (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
    effectiveStress[1] = (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
    effectiveStress[2] = C33 * engineeringStrain[2];

    // calculate the elastic energy density
    const double elasticEnergyDensity = 0.5 * effectiveStress.dot(engineeringStrain);

    double currentEnergyDensity = std::max(elasticEnergyDensity, oldEnergyDensity);

    bool performUpdateAtEnd = false;
    constexpr double residualEnergyDensity = 1.e-12;

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;
            engineeringStress[0] = factor * effectiveStress[0];
            engineeringStress[1] = factor * effectiveStress[1];
            engineeringStress[2] = factor * effectiveStress[2];

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

            // right coefficients are calculated above
            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;

            tangent(0, 0) = factor * C11;
            tangent(1, 0) = factor * C12;
            tangent(2, 0) = 0.;

            tangent(0, 1) = factor * C12;
            tangent(1, 1) = factor * C11;
            tangent(2, 1) = 0.;

            tangent(0, 2) = 0.;
            tangent(1, 2) = 0.;
            tangent(2, 2) = factor * C33;

            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);


            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;

            // plane strain
            engineeringStress3D[0] = factor * (C11 * engineeringStrain[0] + C12 * engineeringStrain[1]);
            engineeringStress3D[1] = factor * (C11 * engineeringStrain[1] + C12 * engineeringStrain[0]);
            engineeringStress3D[2] = factor * C12 * (engineeringStrain[0] + engineeringStrain[1]);
            engineeringStress3D[3] = 0.;
            engineeringStress3D[4] = 0.;
            engineeringStress3D[5] = factor * C33 * engineeringStrain[2];

            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AsEngineeringStrain3D() =
                    engineeringStrain.As3D(mPoissonsRatio, ePlaneState::PLANE_STRAIN);
            break;
        }

        case NuTo::Constitutive::eOutput::ELASTIC_ENERGY_DAMAGED_PART:
        {

            ConstitutiveIOBase& elasticEnergyLoadTerm = *itOutput.second;

            elasticEnergyLoadTerm[0] = currentEnergyDensity;

            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_PHASE_FIELD:
        {

            ConstitutiveIOBase& dStressDPhaseField = *itOutput.second;
            dStressDPhaseField.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            const double degradationFactor = -2 * (1 - damage[0]);


            dStressDPhaseField[0] = degradationFactor * effectiveStress[0];
            dStressDPhaseField[1] = degradationFactor * effectiveStress[1];
            dStressDPhaseField[2] = degradationFactor * effectiveStress[2];


            break;
        }

        case NuTo::Constitutive::eOutput::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
        {
            // the tangent  is equal to the damaged part of the engineering stress for loading and zero for unloading
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);


            if (elasticEnergyDensity == currentEnergyDensity)
            {
                // loading
                tangent[0] = effectiveStress[0];
                tangent[1] = effectiveStress[1];
                tangent[2] = effectiveStress[2];
            }
            else
            {
                // unloading
                tangent.SetZero();
            }

            break;
        }

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        {
            throw Exception(
                    __PRETTY_FUNCTION__,
                    "tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;

        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    // return old/new history variables
    if (performUpdateAtEnd)
    {
        return currentEnergyDensity;
    }
    else
    {
        return oldEnergyDensity;
    }
}


double NuTo::PhaseField::Evaluate2DAnisotropicSpectralDecomposition(const double oldEnergyDensity,
                                                                    const ConstitutiveInputMap& rConstitutiveInput,
                                                                    const ConstitutiveOutputMap& rConstitutiveOutput)
{
    using std::max;
    using std::min;
    using std::pow;
    using std::abs;
    using Eigen::Matrix2d;
    using Eigen::Vector2d;

    constexpr double residualEnergyDensity = 1.e-8;

    const auto& engineeringStrain =
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain2D();
    const auto& damage = *rConstitutiveInput.at(Constitutive::eInput::CRACK_PHASE_FIELD);

    bool performUpdateAtEnd = false;

    Matrix2d strainMatrix;
    strainMatrix(0, 0) = engineeringStrain[0];
    strainMatrix(0, 1) = 0.5 * engineeringStrain[2];
    strainMatrix(1, 0) = 0.5 * engineeringStrain[2];
    strainMatrix(1, 1) = engineeringStrain[1];

    const double traceStrain = strainMatrix.diagonal().sum();

    Eigen::SelfAdjointEigenSolver<Matrix2d> eigenSolver(strainMatrix);

    const double eigenVal00 = eigenSolver.eigenvalues()(0, 0);
    const double eigenVal01 = eigenSolver.eigenvalues()(1, 0);

    const Vector2d eigenVec00 = eigenSolver.eigenvectors().col(0);
    const Vector2d eigenVec01 = eigenSolver.eigenvectors().col(1);

    Matrix2d strainPlus = max(eigenVal00, 0.) * eigenVec00 * eigenVec00.transpose() +
                          max(eigenVal01, 0.) * eigenVec01 * eigenVec01.transpose();

    Matrix2d strainMinus = min(eigenVal00, 0.) * eigenVec00 * eigenVec00.transpose() +
                           min(eigenVal01, 0.) * eigenVec01 * eigenVec01.transpose();

    Matrix2d stressPlus = mLameLambda * max(traceStrain, 0.0) * Matrix2d::Identity() + 2 * mLameMu * strainPlus;

    Matrix2d stressMinus = mLameLambda * min(traceStrain, 0.0) * Matrix2d::Identity() + 2 * mLameMu * strainMinus;

    const double elasticEnergyPlus =
            0.5 * mLameLambda * pow(max(traceStrain, 0.), 2) + mLameMu * (strainPlus * strainPlus).diagonal().sum();

    double currentEnergyDensity = max(elasticEnergyPlus, oldEnergyDensity);

    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor = pow(1. - damage[0], 2) + residualEnergyDensity;

            Matrix2d stressMatrix = factor * stressPlus + stressMinus;

            engineeringStress[0] = stressMatrix(0, 0);
            engineeringStress[1] = stressMatrix(1, 1);
            engineeringStress[2] = stressMatrix(1, 0);

            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsMatrix<3, 3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor = std::pow(1. - damage[0], 2) + residualEnergyDensity;

            CalculateSpectralDecompositionDStressDStrain(tangent, factor, eigenSolver);

            break;
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            throw Exception(__PRETTY_FUNCTION__,
                                     "Visualization of stress not implemented for the anisotropic phase-field model!");
        }

        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            throw Exception(__PRETTY_FUNCTION__,
                                     "Visualization of strain not implemented for the anisotropic phase-field model!");
        }

        case NuTo::Constitutive::eOutput::ELASTIC_ENERGY_DAMAGED_PART:
        {
            ConstitutiveIOBase& elasticEnergyLoadTerm = *itOutput.second;
            elasticEnergyLoadTerm[0] = currentEnergyDensity;

            break;
        }

        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_PHASE_FIELD:
        {

            ConstitutiveIOBase& dStressDPhaseFieldVoigt = *itOutput.second;
            dStressDPhaseFieldVoigt.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            const double factor = -2 * (1 - damage[0]);

            dStressDPhaseFieldVoigt[0] = factor * stressPlus(0, 0);
            dStressDPhaseFieldVoigt[1] = factor * stressPlus(1, 1);
            dStressDPhaseFieldVoigt[2] = factor * stressPlus(1, 0);

            break;
        }

        case NuTo::Constitutive::eOutput::D_ELASTIC_ENERGY_DAMAGED_PART_D_ENGINEERING_STRAIN:
        {

            ConstitutiveIOBase& tangent = *itOutput.second;
            tangent.AssertIsVector<3>(itOutput.first, __PRETTY_FUNCTION__);

            if (elasticEnergyPlus == currentEnergyDensity)
            {
                tangent[0] = stressPlus(0, 0);
                tangent[1] = stressPlus(1, 1);
                tangent[2] = stressPlus(1, 0);
            }
            else
            {
                // unloading
                tangent.SetZero();
            }

            break;
        }

        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
        {
            throw Exception(
                    __PRETTY_FUNCTION__,
                    "tmp_static_data has to be updated without any other outputs, call it separately.");
        }
            continue;

        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {
            performUpdateAtEnd = true;
        }
            continue;

        default:
            continue;
        }
        itOutput.second->SetIsCalculated(true);
    }


    // return old/new history variables
    if (performUpdateAtEnd)
    {
        return currentEnergyDensity;
    }
    else
    {
        return oldEnergyDensity;
    }
}

namespace NuTo
{

template <>
void NuTo::PhaseField::Evaluate<1>(const ConstitutiveInputMap&, const ConstitutiveOutputMap&, Data&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

template <>
void NuTo::PhaseField::Evaluate<2>(const ConstitutiveInputMap& rConstitutiveInput,
                                   const ConstitutiveOutputMap& rConstitutiveOutput, Data& rStaticData)
{
    const auto& planeState =
            *dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get());
    if (planeState.GetPlaneState() != ePlaneState::PLANE_STRAIN)
        throw Exception(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");

    auto itCalculateStaticData = rConstitutiveInput.find(Constitutive::eInput::CALCULATE_STATIC_DATA);
    if (itCalculateStaticData == rConstitutiveInput.end())
        throw Exception(__PRETTY_FUNCTION__,
                                 "You need to specify the way the static data should be calculated (input list).");

    const auto& calculateStaticData =
            dynamic_cast<const ConstitutiveCalculateStaticData&>(*itCalculateStaticData->second);
    int index = calculateStaticData.GetIndexOfPreviousStaticData();

    double oldEnergyDensity = rStaticData.GetData(index);

    double energyDensity = 0.0;

    switch (mEnergyDecomposition)
    {
    case Constitutive::ePhaseFieldEnergyDecomposition::ISOTROPIC:
    {
        energyDensity = Evaluate2DIsotropic(oldEnergyDensity, rConstitutiveInput, rConstitutiveOutput);
        break;
    }
    case Constitutive::ePhaseFieldEnergyDecomposition::ANISOTROPIC_SPECTRAL_DECOMPOSITION:
    {
        energyDensity =
                Evaluate2DAnisotropicSpectralDecomposition(oldEnergyDensity, rConstitutiveInput, rConstitutiveOutput);
        break;
    }
    default:
        throw Exception(__PRETTY_FUNCTION__, "Degradation function type not implemented.");
    }

    // update history variables
    rStaticData.SetData(energyDensity);
}

template <>
void NuTo::PhaseField::Evaluate<3>(const ConstitutiveInputMap&, const ConstitutiveOutputMap&, Data&)
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

} // namespace NuTo

double NuTo::PhaseField::CalculateStaticDataExtrapolationError(ElementBase&, int, const ConstitutiveInputMap&) const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

bool NuTo::PhaseField::CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const
{
    assert(rTimeDerivative == 0 or rTimeDerivative == 1 or rTimeDerivative == 2);
    switch (rTimeDerivative)
    {
    case 0:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::DISPLACEMENTS, Node::eDof::CRACKPHASEFIELD):
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::DISPLACEMENTS):
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
            return true;
        default:
            return false;
        }
        break;
    case 1:
        switch (Node::CombineDofs(rDofRow, rDofCol))
        {
        case Node::CombineDofs(Node::eDof::CRACKPHASEFIELD, Node::eDof::CRACKPHASEFIELD):
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

    switch (rIdentifier)
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
        throw Exception(__PRETTY_FUNCTION__, "Constitutive law does not have the requested variable");
    }
}

void NuTo::PhaseField::SetParameterDouble(Constitutive::eConstitutiveParameter, double)
{
    throw Exception(__PRETTY_FUNCTION__, "Function must not be used");
}

NuTo::Constitutive::eConstitutiveType NuTo::PhaseField::GetType() const
{
    throw Exception(__PRETTY_FUNCTION__, "Not implemented.");
}

void NuTo::PhaseField::CheckParameters() const
{
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
