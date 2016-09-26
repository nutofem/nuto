#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"
#include "nuto/mechanics/interpolationtypes/InterpolationType.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"
#include "nuto/mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "nuto/mechanics/sections/SectionEnum.h"

#define BOOST_TEST_MODULE AdditiveOutputTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

using namespace NuTo;
BOOST_AUTO_TEST_CASE(additive_stresses)
{
    // Create constitutive laws
    AdditiveOutput additiveLaw;

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

    //// mock up an element; ideally, this would not be necessary
    //std::vector<NuTo::NodeBase*> mockNodes;
    //InterpolationType interpolationType(Interpolation::TRIANGLE2D, 2);
    //IntegrationType2D3NGauss1Ip integrationType;
    //interpolationType.UpdateIntegrationType(integrationType);
    //auto element = ContinuumElement<2>(nullptr, mockNodes, ElementData::CONSTITUTIVELAWIP, IpData::NOIPDATA, &interpolationType);
    //auto section = SectionPlane(Section::PLANE_STRESS);
    //element.SetSection(&section);

    // evaluate the additive input law
    additiveLaw.Evaluate2D(inputMap, outputMap, nullptr);

    // compare to expected results
    const auto& stress = *static_cast<EngineeringStress<2>*>(outputMap.at(Constitutive::eOutput::ENGINEERING_STRESS).get());
    Eigen::Matrix<double, 3, 1> expectedStress;
    expectedStress << 2.0, 2.0, 0.0;
    BOOST_CHECK_SMALL((stress - expectedStress).norm(), 1e-16);
}
