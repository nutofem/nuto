#include "nuto/mechanics/constitutive/laws/AdditiveInputExplicit.h"
#include "nuto/mechanics/constitutive/laws/LinearElasticEngineeringStress.h"
#include "nuto/mechanics/constitutive/laws/ThermalStrains.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/sections/SectionPlane.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D3NGauss1Ip.h"

#define BOOST_TEST_MODULE AdditiveInputTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

using namespace NuTo;
BOOST_AUTO_TEST_CASE(additive_strains)
{
    // Create constitutive laws
    AdditiveInputExplicit additiveLaw;

    LinearElasticEngineeringStress linElastic;
    linElastic.SetParameterDouble(Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 1.0);
    linElastic.SetParameterDouble(Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);

    ThermalStrains thermalStrains;
    thermalStrains.SetParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_EXPANSION_COEFFICIENT, 0.5);

    additiveLaw.AddConstitutiveLaw(linElastic);
    additiveLaw.AddConstitutiveLaw(thermalStrains, Constitutive::Input::ENGINEERING_STRAIN);

    // Create input data
    ConstitutiveInputMap inputMap;
    inputMap[Constitutive::Input::TEMPERATURE] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Input::TEMPERATURE);
    inputMap[Constitutive::Input::ENGINEERING_STRAIN] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Input::ENGINEERING_STRAIN);
    inputMap[Constitutive::Input::PLANE_STATE] = ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Input::PLANE_STATE);
    (*static_cast<ConstitutiveScalar*>(inputMap.at(Constitutive::Input::TEMPERATURE).get()))[0] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::Input::ENGINEERING_STRAIN).get()))[0] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::Input::ENGINEERING_STRAIN).get()))[1] = 1.0;
    (*static_cast<EngineeringStrain<2>*>(inputMap.at(Constitutive::Input::ENGINEERING_STRAIN).get()))[2] = 0.0;

    // ask for stress as output
    ConstitutiveOutputMap outputMap;
    outputMap[Constitutive::Output::ENGINEERING_STRESS] = 
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Output::ENGINEERING_STRESS);
    outputMap[Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN] = 
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Output::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN);
    outputMap[Constitutive::Output::THERMAL_STRAIN] =
        ConstitutiveIOBase::makeConstitutiveIO<2>(Constitutive::Output::THERMAL_STRAIN);

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
    const auto& stress = *static_cast<EngineeringStress<2>*>(outputMap.at(Constitutive::Output::ENGINEERING_STRESS).get());
    Eigen::Matrix<double, 3, 1> expectedStress;
    expectedStress << 0.5, 0.5, 0.0;
    BOOST_CHECK_SMALL((stress - expectedStress).norm(), 1e-16);
}
