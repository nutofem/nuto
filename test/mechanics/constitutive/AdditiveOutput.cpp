#include "BoostUnitTest.h"

#include "mechanics/constitutive/laws/AdditiveOutput.h"
#include "mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "mechanics/constitutive/laws/ThermalStrains.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"


using namespace NuTo;
BOOST_AUTO_TEST_CASE(additive_stresses)
{
    int NumTimeDerivatives = 2;
    // Create constitutive laws
    AdditiveOutput additiveLaw(NumTimeDerivatives);

    // two "springs" in parallel; their stiffnesses add up
    // i.e. unit strain should result in stresses of two
    LinearElasticEngineeringStress linElasticOne;
    linElasticOne.SetParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
    linElasticOne.SetParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    LinearElasticEngineeringStress linElasticTwo;
    linElasticTwo.SetParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
    linElasticTwo.SetParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    additiveLaw.AddConstitutiveLaw(linElasticOne);
    additiveLaw.AddConstitutiveLaw(linElasticTwo);

    // Create input data
    ConstitutiveInputMap inputMap;

    inputMap[Constitutive::eInput::ENGINEERING_STRAIN] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eInput::ENGINEERING_STRAIN);
    inputMap[Constitutive::eInput::PLANE_STATE] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eInput::PLANE_STATE);
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[0] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[1] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::eInput::ENGINEERING_STRAIN).get()))[2] = 0.0;

    // ask for stress as output
    ConstitutiveOutputMap outputMap;
    outputMap[Constitutive::eOutput::ENGINEERING_STRESS] =
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::eOutput::ENGINEERING_STRESS);

    // evaluate the additive input law
    auto additiveIP = additiveLaw.CreateIPLaw();
    additiveIP->Evaluate<2>(inputMap, outputMap);

    // compare to expected results
    const auto& stress = *static_cast<EngineeringStress<2>*>(outputMap.at(Constitutive::eOutput::ENGINEERING_STRESS).get());
    Eigen::Matrix<double, 3, 1> expectedStress;
    expectedStress << 2.0, 2.0, 0.0;
    BOOST_CHECK_SMALL((stress - expectedStress).norm(), 1e-16);
}
