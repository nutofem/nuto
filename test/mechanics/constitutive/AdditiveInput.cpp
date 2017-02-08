#include "BoostUnitTest.h"

#include "mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/sections/SectionEnum.h"

using namespace NuTo;

BOOST_AUTO_TEST_CASE(additive_strains)
{
    int NumTimeDerivatives = 2;
    // Create constitutive laws
    AdditiveInputExplicit additiveLaw(NumTimeDerivatives);

    LinearElasticEngineeringStress linElastic;
    linElastic.SetParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
    linElastic.SetParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    ThermalStrains thermalStrains;
    thermalStrains.SetParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 0.5);

    additiveLaw.AddConstitutiveLaw(linElastic);
    additiveLaw.AddConstitutiveLaw(thermalStrains, Constitutive::eInput::ENGINEERING_STRAIN);

    // Create input data
    ConstitutiveInputMap inputMap;
    inputMap[Constitutive::eInput::TEMPERATURE] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eInput::TEMPERATURE);
    inputMap[Constitutive::eInput::ENGINEERING_STRAIN] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eInput::ENGINEERING_STRAIN);
    inputMap[Constitutive::eInput::PLANE_STATE] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eInput::PLANE_STATE);
    (*static_cast<ConstitutiveScalar*>(inputMap.at(Constitutive::eInput::TEMPERATURE).get()))[0] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[0] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[1] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[2] = 0.0;

    // ask for stress as output
    ConstitutiveOutputMap outputMap;
    outputMap[Constitutive::eOutput::ENGINEERING_STRESS] =
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eOutput::ENGINEERING_STRESS);
    outputMap[Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] =
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);
    outputMap[Constitutive::eOutput::THERMAL_STRAIN] =
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eOutput::THERMAL_STRAIN);

    // evaluate the additive input law
    auto additiveIP = additiveLaw.CreateIPLaw();
    additiveIP->Evaluate<2>(inputMap, outputMap);

    // compare to expected results
    const auto& stress = *static_cast<EngineeringStress<2>*>(outputMap.at(Constitutive::eOutput::ENGINEERING_STRESS).get());
    Eigen::Matrix<double, 3, 1> expectedStress;
    expectedStress << 0.5, 0.5, 0.0;
    BOOST_CHECK_SMALL((stress - expectedStress).norm(), 1e-16);
}
