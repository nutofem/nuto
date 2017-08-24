#include "mechanics/constitutive/laws/Creep.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/laws/EngineeringStressHelper.h"
#include "mechanics/nodes/NodeEnum.h"


using namespace NuTo;
using namespace NuTo::Constitutive;

Creep::Creep()
    : ConstitutiveBase()
{
}


template <int TDim>
void Creep::Evaluate(const ConstitutiveInputMap& rConstitutiveInput, const ConstitutiveOutputMap& rConstitutiveOutput,
                     Creep::Data& rStaticData)
{
    assert(mKC_E.rows() > 0 || mE > 0.0);
    assert(mKC_E.cols() == 1);
    assert(mKC_E.rows() == mKC_T.rows());
    assert(mKC_E.cols() == mKC_T.cols());

    const unsigned int numKelvinUnits = mKC_E.rows();

    constexpr int VoigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);

    // get static data
    auto& staticData = rStaticData.GetData();

    // get and set time variables
    const auto& itTime = rConstitutiveInput.find(Constitutive::eInput::TIME);
    if (itTime != rConstitutiveInput.end())
        staticData.SetCurrentTime((*itTime->second)[0]);
    double delta_t = staticData.GetCurrentTime() - staticData.GetPreviousTime();

    // initialize history data
    if (staticData.mHistoryData.rows() != VoigtDim || staticData.mHistoryData.cols() != mKC_E.rows())
    {
        staticData.mHistoryData = Eigen::MatrixXd::Zero(VoigtDim, mKC_E.rows());
        staticData.mHistoryStrain = Eigen::VectorXd::Zero(VoigtDim);
        staticData.mHistoryStress = Eigen::VectorXd::Zero(VoigtDim);
        staticData.mDeltaStress = Eigen::VectorXd::Zero(VoigtDim);
        staticData.mDeltaStrain = Eigen::VectorXd::Zero(VoigtDim);
        staticData.mDeltaCreepStrain = Eigen::VectorXd::Zero(VoigtDim);
    }
    assert(staticData.mHistoryData.rows() == VoigtDim);
    assert(staticData.mHistoryData.cols() == numKelvinUnits);
    assert(staticData.mHistoryStrain.rows() == VoigtDim);
    assert(staticData.mHistoryStrain.cols() == 1);
    assert(staticData.mHistoryStress.rows() == VoigtDim);
    assert(staticData.mHistoryStress.cols() == 1);
    assert(staticData.mDeltaStress.rows() == VoigtDim);
    assert(staticData.mDeltaStress.cols() == 1);
    assert(staticData.mDeltaStrain.rows() == VoigtDim);
    assert(staticData.mDeltaStrain.cols() == 1);
    assert(staticData.mDeltaCreepStrain.rows() == VoigtDim);
    assert(staticData.mDeltaCreepStrain.cols() == 1);


    for (auto& itOutput : rConstitutiveOutput)
    {
        switch (itOutput.first)
        {
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS:
        {
            // Get and check necessary references
            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)
                    ->AssertIsVector<VoigtDim>(itOutput.first, __PRETTY_FUNCTION__);
            Eigen::Matrix<double, VoigtDim, 1>& engineeringStrain =
                    static_cast<ConstitutiveVector<VoigtDim>*>(
                            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN).get())
                            ->AsVector();
            ConstitutiveIOBase& engineeringStress = *itOutput.second;
            engineeringStress.AssertIsVector<VoigtDim>(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            staticData.mDeltaCreepStrain = Eigen::VectorXd::Zero(VoigtDim);
            staticData.mDeltaStrain = 1. * (engineeringStrain - staticData.mHistoryStrain);

            for (unsigned int i = 0; i < VoigtDim; ++i)
            {
                //                engineeringStress[i] = mKC_E(0) * engineeringStrain[i]; // still linear elastic
                for (unsigned int j = 0; j < numKelvinUnits; ++j)
                {
                    staticData.mDeltaCreepStrain[i] +=
                            (1 - ExponentialAlgorithmCalculateBeta(j, delta_t)) * staticData.mHistoryData(i, j);
                }
            }
            staticData.mDeltaStress = ExponentialAlgorithmCalculateStiffnessMatrix<TDim>(delta_t, rConstitutiveInput) *
                                      (staticData.mDeltaStrain - staticData.mDeltaCreepStrain);
            for (unsigned int i = 0; i < VoigtDim; ++i)
            {
                engineeringStress[i] = staticData.mHistoryStress[i] + staticData.mDeltaStress[i];
            }
            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRESS_VISUALIZE:
        {
            // Get and check necessary references
            ConstitutiveIOBase& engineeringStress3D = *itOutput.second;
            engineeringStress3D.AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);

            rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)
                    ->AssertIsVector<VoigtDim>(itOutput.first, __PRETTY_FUNCTION__);

            // Calculation
            engineeringStress3D.SetZero();
            for (unsigned int i = 0; i < VoigtDim; ++i)
            {
                engineeringStress3D[i] = staticData.mHistoryStress[i];
            }
            break;
        }
        case NuTo::Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN:
        {
            // Get and check necessary references

            dynamic_cast<Eigen::Matrix<double, VoigtDim, VoigtDim>&>(*itOutput.second) =
                    ExponentialAlgorithmCalculateStiffnessMatrix<TDim>(delta_t, rConstitutiveInput);


            break;
        }
        case NuTo::Constitutive::eOutput::ENGINEERING_STRAIN_VISUALIZE:
        {
            itOutput.second->AssertIsVector<6>(itOutput.first, __PRETTY_FUNCTION__);
            switch (TDim)
            {
            case 1:
                itOutput.second->AsEngineeringStrain3D() =
                        rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)
                                ->AsEngineeringStrain1D()
                                .As3D(mNu);
                break;
            case 2:
                itOutput.second->AsEngineeringStrain3D() =
                        rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)
                                ->AsEngineeringStrain2D()
                                .As3D(mNu);
                break;
            case 3:
                itOutput.second->AsEngineeringStrain3D() =
                        rConstitutiveInput.at(Constitutive::eInput::ENGINEERING_STRAIN)->AsEngineeringStrain3D();
                break;
            }
            break;
        }
        case NuTo::Constitutive::eOutput::UPDATE_STATIC_DATA:
        {

            const double KelvinChainStiffness = ExponentialAlgorithmCalculateChainStiffness(delta_t);

            for (unsigned int j = 0; j < numKelvinUnits; ++j)
            {
                double beta = ExponentialAlgorithmCalculateBeta(j, delta_t);
                double lambda = ExponentialAlgorithmCalculateLambda(j, delta_t);
                for (unsigned int i = 0; i < VoigtDim; ++i)
                {
                    staticData.mHistoryData(i, j) =
                            lambda * KelvinChainStiffness / mKC_E(j) *
                                    (staticData.mDeltaStrain(i) - staticData.mDeltaCreepStrain(i)) +
                            beta * staticData.mHistoryData(i, j);
                }
            }

            staticData.ProceedToNextTimestep();
            staticData.mHistoryStress += staticData.mDeltaStress;
            staticData.mHistoryStrain += staticData.mDeltaStrain;
            staticData.mDeltaStress = Eigen::VectorXd::Zero(VoigtDim);
            staticData.mDeltaStrain = Eigen::VectorXd::Zero(VoigtDim);
            staticData.mDeltaCreepStrain = Eigen::VectorXd::Zero(VoigtDim);
        }
            continue;
        case NuTo::Constitutive::eOutput::EXTRAPOLATION_ERROR:
            break;
        case NuTo::Constitutive::eOutput::UPDATE_TMP_STATIC_DATA:
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


ConstitutiveInputMap NuTo::Creep::GetConstitutiveInputs(const ConstitutiveOutputMap& rConstitutiveOutput) const
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
        }
    }

    return constitutiveInputMap;
}


bool NuTo::Creep::CheckDofCombinationComputable(Node::eDof rDofRow, Node::eDof rDofCol, int rTimeDerivative) const
{
    if (rTimeDerivative < 1 && rDofCol == Node::eDof::DISPLACEMENTS && rDofRow == Node::eDof::DISPLACEMENTS)
        return true;
    else
        return false;
}

Constitutive::eConstitutiveType NuTo::Creep::GetType() const
{
    return Constitutive::eConstitutiveType::CREEP;
}

bool NuTo::Creep::HaveTmpStaticData() const
{
    return false;
}

void NuTo::Creep::CheckParameters() const
{
}

void Creep::SetParameterDouble(Constitutive::eConstitutiveParameter identifier, double value)
{
    switch (identifier)
    {
    case eConstitutiveParameter::YOUNGS_MODULUS:
        mE = value;
        break;
    case eConstitutiveParameter::POISSONS_RATIO:
        mNu = value;
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(identifier));
    }
}

void Creep::SetParameterFullVectorDouble(Constitutive::eConstitutiveParameter identifier, Eigen::VectorXd value)
{
    switch (identifier)
    {
    case eConstitutiveParameter::KELVIN_CHAIN_DAMPING:
        mKC_D = value;
        break;
    case eConstitutiveParameter::KELVIN_CHAIN_RETARDATIONTIME:
        mKC_T = value;
        break;
    case eConstitutiveParameter::KELVIN_CHAIN_STIFFNESS:
        mKC_E = value;
        break;
    default:
        throw Exception(__PRETTY_FUNCTION__, std::string("Constitutive law does not have the parameter ") +
                                                     Constitutive::ConstitutiveParameterToString(identifier));
    }
}

template <int TDim>
Eigen::MatrixXd
Creep::ExponentialAlgorithmCalculateStiffnessMatrix(double delta_t,
                                                    const ConstitutiveInputMap& rConstitutiveInput) const
{
    constexpr int voigtDim = ConstitutiveIOBase::GetVoigtDim(TDim);
    Eigen::MatrixXd stiffnessMat = Eigen::MatrixXd::Zero(voigtDim, voigtDim);

    switch (TDim)
    {
    case 1:
        stiffnessMat(0, 0) = ExponentialAlgorithmCalculateChainStiffness(delta_t);
        break;


    case 2:
    {
        double C11, C12, C33;
        switch ((*dynamic_cast<ConstitutivePlaneState*>(rConstitutiveInput.at(Constitutive::eInput::PLANE_STATE).get()))
                        .GetPlaneState())
        {
        case ePlaneState::PLANE_STRAIN:
            std::tie(C11, C12, C33) = NuTo::EngineeringStressHelper::CalculateCoefficients3D(
                    ExponentialAlgorithmCalculateChainStiffness(delta_t), mNu);
            break;
        case ePlaneState::PLANE_STRESS:
            std::tie(C11, C12, C33) = EngineeringStressHelper::CalculateCoefficients2DPlaneStress(
                    ExponentialAlgorithmCalculateChainStiffness(delta_t), mNu);
            break;
        default:
            throw Exception(__PRETTY_FUNCTION__, "Invalid type of 2D section behavior found.");
        }

        stiffnessMat(0, 0) = C11;
        stiffnessMat(1, 0) = C12;
        stiffnessMat(0, 1) = C12;
        stiffnessMat(1, 1) = C11;
        stiffnessMat(2, 2) = C33;
    }
    break;


    case 3:
    {
        double C11, C12, C44;
        std::tie(C11, C12, C44) = EngineeringStressHelper::CalculateCoefficients3D(
                ExponentialAlgorithmCalculateChainStiffness(delta_t), mNu);

        // C11 diagonal:
        stiffnessMat(0, 0) = C11;
        stiffnessMat(1, 1) = C11;
        stiffnessMat(2, 2) = C11;
        // C12 off diagonals:
        stiffnessMat(0, 1) = C12;
        stiffnessMat(0, 2) = C12;
        stiffnessMat(1, 0) = C12;
        stiffnessMat(1, 2) = C12;
        stiffnessMat(2, 0) = C12;
        stiffnessMat(2, 1) = C12;
        // C44 diagonal:
        stiffnessMat(3, 3) = C44;
        stiffnessMat(4, 4) = C44;
        stiffnessMat(5, 5) = C44;
    }
    break;


    default:
        throw Exception(__PRETTY_FUNCTION__, "Invalid dimension!");
    }


    return stiffnessMat;
}

double Creep::ExponentialAlgorithmCalculateBeta(unsigned int index, double delta_t) const
{
    assert(index < mKC_T.rows());
    return std::exp(-delta_t / mKC_T[index]);
}

double Creep::ExponentialAlgorithmCalculateLambda(unsigned int index, double delta_t) const
{
    assert(index < mKC_T.rows());
    if (delta_t > 0.0) // <--- should be an assert, but that conflicts with Newmarks initial state calculation
        return mKC_T[index] / delta_t * (1 - ExponentialAlgorithmCalculateBeta(index, delta_t));
    else
        return 0.;
}

double Creep::ExponentialAlgorithmCalculateChainStiffness(double delta_t) const
{
    assert(mKC_E.rows() > 0 || mE > 0.0);
    assert(mKC_E.cols() == 1);
    double KelvinChainStiffnessInverted = 0.0;
    if (mE > 0.)
        KelvinChainStiffnessInverted += 1. / mE;
    for (unsigned int i = 0; i < mKC_E.rows(); ++i)
    {
        assert(mKC_E[i] > 0.0);
        KelvinChainStiffnessInverted += (1. - ExponentialAlgorithmCalculateLambda(i, delta_t)) / mKC_E[i];
    }
    return 1. / KelvinChainStiffnessInverted;
}
