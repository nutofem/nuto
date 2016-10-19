
#include "nuto/mechanics/constitutive/laws/LocalDamageModel.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/staticData/Leaf.h"


#include <iostream>

#define BOOST_TEST_MODULE LocalDamageModelTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}


using std::cout;
using std::endl;
using NuTo::Constitutive::eInput;
using NuTo::Constitutive::eOutput;
using NuTo::Constitutive::eConstitutiveParameter;

double EvaluateLocalDamageModelModel(const Eigen::Vector3d& rStrain, NuTo::LocalDamageModel* rLocalDamageModel)
{
    constexpr int VoigtDim = NuTo::ConstitutiveIOBase::GetVoigtDim(2);

    NuTo::ConstitutiveInputMap myConstitutiveInputMap;
    myConstitutiveInputMap[eInput::CALCULATE_STATIC_DATA] = std::make_unique<NuTo::ConstitutiveCalculateStaticData>(NuTo::eCalculateStaticData::EULER_BACKWARD);

    myConstitutiveInputMap[eInput::ENGINEERING_STRAIN]    = NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eInput::ENGINEERING_STRAIN);
    auto& engineeringStrain= *static_cast<NuTo::EngineeringStrain<2>*>(myConstitutiveInputMap.at(eInput::ENGINEERING_STRAIN).get());
    engineeringStrain.AsVector() = rStrain;

    cout << "Inputs:"                         << endl;
    cout << "engineeringStrain.AsVector()"    << endl << engineeringStrain.AsVector()   << endl;

    NuTo::ConstitutiveOutputMap myConstitutiveOutputMap;
    myConstitutiveOutputMap[eOutput::ENGINEERING_STRESS] = NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eOutput::ENGINEERING_STRESS);
    auto& engineeringStress = *static_cast<NuTo::EngineeringStress<2>*>(myConstitutiveOutputMap.at(eOutput::ENGINEERING_STRESS).get());
    engineeringStress.SetZero();

    myConstitutiveOutputMap[eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] = NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);
    auto& tangentStressStrain = *static_cast<NuTo::ConstitutiveMatrix<VoigtDim, VoigtDim>*>(myConstitutiveOutputMap.at(eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN).get());
    tangentStressStrain.SetZero();

    NuTo::Constitutive::StaticData::Leaf<double> staticData;
    staticData.SetData(0.0);

    rLocalDamageModel->Evaluate2D(NuTo::eSectionType::PLANE_STRAIN, staticData, myConstitutiveInputMap, myConstitutiveOutputMap);

    cout << "Outputs:"                       << endl;
    cout << "engineeringStress.AsVector()"   << endl << engineeringStress.AsVector()   << endl;
    cout << "tangentStressStrain.AsMatrix()" << endl << tangentStressStrain.AsMatrix() << endl;

    Eigen::Vector3d residual = tangentStressStrain.AsMatrix() * engineeringStrain - engineeringStress;

    return residual.norm();
}

BOOST_AUTO_TEST_CASE(check_d_stress_d_strain)
{
    constexpr double tolerance              = 1.e-10;
    constexpr double youngsModulus          = 4e4;
    constexpr double poissonsRatio          = 0.2;
    constexpr double tensileStrength        = 3;
    constexpr double compressiveStrength    = 30;
    constexpr double fractureEnergy         = 0.01;


    NuTo::LocalDamageModel* localDamageModel = new NuTo::LocalDamageModel();
    localDamageModel->SetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS,       youngsModulus);
    localDamageModel->SetParameterDouble(eConstitutiveParameter::POISSONS_RATIO,       poissonsRatio);
    localDamageModel->SetParameterDouble(eConstitutiveParameter::TENSILE_STRENGTH,     tensileStrength);
    localDamageModel->SetParameterDouble(eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    localDamageModel->SetParameterDouble(eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy);

    Eigen::Vector3d strain;
    strain[0] = 0.;
    strain[1] = 0.;
    strain[2] = 0.;


    double residualNorm = EvaluateLocalDamageModelModel(strain, localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.e-5;
    strain[1] = 0.;
    strain[2] = 0.;

    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = -1.e-5;
    strain[1] = 0.;
    strain[2] = 0.;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.e-5;
    strain[1] = 1.e-5;
    strain[2] = 0.;

    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 2.e-5;
    strain[1] = 1.e-5;
    strain[2] = 0.;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 2.e-5;
    strain[1] = -1.e-5;
    strain[2] = 0.;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 0.;
    strain[1] = 0.;
    strain[2] = 2.e-5;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;
    strain[0] = 1.e-5;
    strain[1] = 1.e-5;
    strain[2] = 2.e-5;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);

    cout << "**********************************************" << endl;

    strain[0] = 1.e-5;
    strain[1] = -3.e-5;
    strain[2] = 2.e-5;


    residualNorm = EvaluateLocalDamageModelModel(strain,  localDamageModel);
    BOOST_REQUIRE_SMALL(residualNorm,tolerance);
}


