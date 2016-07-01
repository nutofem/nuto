#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/constitutive/laws/PhaseField.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/elements/EvaluateDataContinuum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataHistoryVariableScalar.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <exception>
#include <cassert>

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
void EvaluatePhaseFieldModel(const Eigen::Vector3d& rStrain, const double rCrackPhaseField, NuTo::PhaseField* rPhaseField)
{


    NuTo::EvaluateDataContinuum<2> data;
    data.mEngineeringStress.AsVector().setZero();

    NuTo::ConstitutiveInputMap myConstitutiveInputMap;
    myConstitutiveInputMap[NuTo::Constitutive::Input::CALCULATE_STATIC_DATA] = std::make_unique<NuTo::ConstitutiveCalculateStaticData>(NuTo::CalculateStaticData::EULER_BACKWARD);
    myConstitutiveInputMap[NuTo::Constitutive::Input::ENGINEERING_STRAIN]    = NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::Input::ENGINEERING_STRAIN);
    auto& engineeringStrain= *static_cast<NuTo::EngineeringStrain<2>*>(myConstitutiveInputMap.at(NuTo::Constitutive::Input::ENGINEERING_STRAIN).get());
    engineeringStrain.AsVector() = rStrain;
    myConstitutiveInputMap[NuTo::Constitutive::Input::CRACK_PHASE_FIELD]     = NuTo::ConstitutiveIOBase::makeConstitutiveIO<2>(NuTo::Constitutive::Input::CRACK_PHASE_FIELD);
    auto& damage = *static_cast<NuTo::ConstitutiveScalar*>(myConstitutiveInputMap.at(NuTo::Constitutive::Input::CRACK_PHASE_FIELD).get());
    damage[0] = rCrackPhaseField;

    cout << "Inputs:" << endl;
    cout << "engineeringStrain.AsVector()"    << endl << engineeringStrain.AsVector()   << endl;
    cout << "damage.AsScalar()"               << endl << damage.AsScalar()              << endl;


    NuTo::ConstitutiveStaticDataHistoryVariableScalar staticData;
    staticData.SetHistoryVariable(0.0);

    NuTo::ConstitutiveOutputMap myConstitutiveOutputMap;
    myConstitutiveOutputMap[NuTo::Constitutive::Output::ENGINEERING_STRESS]                            = &data.mEngineeringStress;
    myConstitutiveOutputMap[NuTo::Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN]     = &data.mTangentStressStrain;

    rPhaseField->Evaluate2DAnisotropicSpectralDecomposition(staticData, myConstitutiveInputMap, myConstitutiveOutputMap);

    cout << "Outputs:" << endl;
    cout << "data.mEngineeringStress.AsVector()"    << endl << data.mEngineeringStress.AsVector()   << endl;
    cout << "data.mTangentStressStrain.AsMatrix()"    << endl << data.mTangentStressStrain.AsMatrix()   << endl;

    Eigen::Vector3d residual = data.mTangentStressStrain.AsMatrix() * engineeringStrain - data.mEngineeringStress;
    cout << "Residual"    << endl << residual   << endl;

    assert(residual.norm() < 1e-8);
}


int main()
{
    try
    {

        constexpr double youngsModulus = 200;
        constexpr double poissonsRatio = 0.3;
        constexpr double lengthScale = 1.0;
        constexpr double fractureEnergy = 2.7;
        constexpr double artificialViscosity = 0.1;
        constexpr   ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ANISOTROPIC_SPECTRAL_DECOMPOSITION;


        cout << "**********************************************" << endl;
        cout << "**  material                                **" << endl;
        cout << "**********************************************" << endl;

        NuTo::PhaseField* phaseField = new NuTo::PhaseField(youngsModulus,
                                                            poissonsRatio,
                                                            lengthScale,
                                                            fractureEnergy,
                                                            artificialViscosity,
                                                            energyDecomposition);

        cout << "**********************************************" << endl;
        Eigen::Vector3d strain;
        strain[0] = 0.;
        strain[1] = 0.;
        strain[2] = 0.;

        double crackPhaseField = 0.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 0.;
        strain[1] = 0.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);


        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = 0.;
        strain[2] = 0.;

        crackPhaseField = 0.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);
        cout << "**********************************************" << endl;
        strain[0] = -1.;
        strain[1] = 0.;
        strain[2] = 0.;

        crackPhaseField = 0.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = -1.;
        strain[1] = 0.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = 0.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = 1.;
        strain[2] = 0.;

        crackPhaseField = 0.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = 1.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);


        cout << "**********************************************" << endl;
        strain[0] = 2.;
        strain[1] = 1.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 2.;
        strain[1] = -1.;
        strain[2] = 0.;

        crackPhaseField = 1.;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);


        cout << "**********************************************" << endl;
        strain[0] = 0.;
        strain[1] = 0.;
        strain[2] = 2.;

        crackPhaseField = 1.0;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);


        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = 1.;
        strain[2] = 2.;

        crackPhaseField = 1.0;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

        cout << "**********************************************" << endl;
        strain[0] = 1.;
        strain[1] = -3.;
        strain[2] = 2.;

        crackPhaseField = 1.0;

        EvaluatePhaseFieldModel(strain, crackPhaseField, phaseField);

    } catch (...)
    {
        cout << "Test failed" << endl;
        return EXIT_FAILURE;
    }

    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
