#include "BoostUnitTest.h"
#include "ConstitutiveTangentTester.h"

#include "mechanics/constitutive/laws/LocalDamageModel.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"

#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"

#include <iostream>

using std::cout;
using std::endl;
using NuTo::Constitutive::eInput;
using NuTo::Constitutive::eOutput;
using NuTo::Constitutive::eConstitutiveParameter;

template <int TDim>
void EvaluateLocalDamageModelModel(NuTo::EngineeringStrain<TDim> rStrain, NuTo::LocalDamageModel& rLocalDamageModel,
                                   double rKappa)
{
    auto iplaw = rLocalDamageModel.CreateIPLaw();
    iplaw->GetData<NuTo::LocalDamageModel>().SetData(rKappa);
    NuTo::Test::ConstitutiveTangentTester<TDim> tester(*iplaw.get(), 1.e-8, 2.e-5);

    NuTo::ConstitutiveInputMap input;
    input.Add<TDim>(eInput::ENGINEERING_STRAIN);
    input.Add<TDim>(eInput::CALCULATE_STATIC_DATA);
    input.Add<TDim>(eInput::PLANE_STATE);

    input[eInput::ENGINEERING_STRAIN]->AsEngineeringStrain<TDim>() = rStrain;
    dynamic_cast<NuTo::ConstitutiveCalculateStaticData&>(*input[eInput::CALCULATE_STATIC_DATA])
            .SetCalculateStaticData(NuTo::eCalculateStaticData::EULER_BACKWARD);

    auto& planeState = dynamic_cast<NuTo::ConstitutivePlaneState&>(*input[eInput::PLANE_STATE]);
    planeState.SetPlaneState(NuTo::ePlaneState::PLANE_STRESS);


    cout << "Input strain: " << rStrain.CopyToEigenMatrix().transpose() << endl;


    BOOST_CHECK(tester.CheckTangent(input, eInput::ENGINEERING_STRAIN, eOutput::ENGINEERING_STRESS,
                                    eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));

    if (TDim == 2)
    {
        planeState.SetPlaneState(NuTo::ePlaneState::PLANE_STRAIN);

        BOOST_CHECK(tester.CheckTangent(input, eInput::ENGINEERING_STRAIN, eOutput::ENGINEERING_STRESS,
                                        eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));
    }
}

auto GetLocalDamageModel()
{
    constexpr double youngsModulus = 4e4;
    constexpr double poissonsRatio = 0.2;
    constexpr double tensileStrength = 3;
    constexpr double compressiveStrength = 30;
    constexpr double fractureEnergy = 0.01;

    NuTo::LocalDamageModel localDamageModel;
    localDamageModel.SetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    localDamageModel.SetParameterDouble(eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    localDamageModel.SetParameterDouble(eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    localDamageModel.SetParameterDouble(eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    localDamageModel.SetDamageLaw(NuTo::Constitutive::DamageLawExponential::Create(tensileStrength / youngsModulus,
                                                                                   tensileStrength / fractureEnergy));
    return localDamageModel;
}

BOOST_AUTO_TEST_CASE(check_d_stress_d_strain2D)
{
    NuTo::LocalDamageModel localDamageModel = GetLocalDamageModel();
    double kappa_0 = 3. / 4e4;
    double kappa = kappa_0 / 3.;
    EvaluateLocalDamageModelModel<2>({0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({1.e-5, 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({-1.e-5, 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({1.e-5, 1.e-5, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({2.e-5, 1.e-5, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({2.e-5, -1.e-5, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({0, 0, 2.e-5}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({1.e-5, 1.e-5, 2.e-5}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({1.e-5, -2.e-5, 2.e-5}, localDamageModel, kappa);

    // some test in damaged loading
    kappa = 2 * kappa_0;
    double eps = 1.e-7; // small load increment = damaged loading
    EvaluateLocalDamageModelModel<2>({kappa + eps, 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, eps, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, 0., eps}, localDamageModel, kappa);

    EvaluateLocalDamageModelModel<2>({kappa + eps, +eps, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, eps, eps}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa + eps, 0., eps}, localDamageModel, kappa);


    // decrement = elastic unloading
    EvaluateLocalDamageModelModel<2>({kappa - eps, 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, -eps, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, 0., -eps}, localDamageModel, kappa);

    EvaluateLocalDamageModelModel<2>({kappa - eps, -eps, 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa, -eps, -eps}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<2>({kappa - eps, 0., -eps}, localDamageModel, kappa);
}


BOOST_AUTO_TEST_CASE(check_d_stress_d_strain3D)
{
    NuTo::LocalDamageModel localDamageModel = GetLocalDamageModel();
    double kappa_0 = 3. / 4e4;
    double kappa = kappa_0 / 3.;
    EvaluateLocalDamageModelModel<3>({0., 0., 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({1.e-5, 0., 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({-1.e-5, 0., 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({1.e-5, 1.e-5, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({2.e-5, 1.e-5, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({2.e-5, -1.e-5, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({0, 0, 2.e-5, 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({1.e-5, 1.e-5, 2.e-5, 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({1.e-5, -2.e-5, 2.e-5, 0., 0., 0.}, localDamageModel, kappa);

    // some test in damaged loading
    kappa = 2 * kappa_0;
    double eps = 1.e-6; // small load increment = damaged loading
    EvaluateLocalDamageModelModel<3>({kappa + eps, 0., 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, eps, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, 0., eps, 0., 0., 0.}, localDamageModel, kappa);

    EvaluateLocalDamageModelModel<3>({kappa + eps, +eps, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, eps, eps, 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa + eps, 0., eps, 0., 0., 0.}, localDamageModel, kappa);


    // decrement = elastic unloading
    EvaluateLocalDamageModelModel<3>({kappa - eps, 0., 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, -eps, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, 0., -eps, 0., 0., 0.}, localDamageModel, kappa);

    EvaluateLocalDamageModelModel<3>({kappa - eps, -eps, 0., 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa, -eps, -eps, 0., 0., 0.}, localDamageModel, kappa);
    EvaluateLocalDamageModelModel<3>({kappa - eps, 0., -eps, 0., 0., 0.}, localDamageModel, kappa);
}

BOOST_AUTO_TEST_CASE(LocalDamageModelVisualize3D)
{
    NuTo::LocalDamageModel localDamageModel;
    localDamageModel.SetDamageLaw(NuTo::Constitutive::DamageLawExponential::Create(1, 2));
    auto law = localDamageModel.CreateIPLaw();
    NuTo::ConstitutiveInputMap input;
    input.Add<3>(eInput::ENGINEERING_STRAIN);
    input.Add<3>(eInput::CALCULATE_STATIC_DATA);

    NuTo::ConstitutiveOutputMap output;
    output.Add<3>(eOutput::ENGINEERING_STRAIN_VISUALIZE);
    output.Add<3>(eOutput::ENGINEERING_STRESS_VISUALIZE);
    output.Add<3>(eOutput::DAMAGE);
    output.Add<3>(eOutput::LOCAL_EQ_STRAIN);
    law->Evaluate<3>(input, output);
}
