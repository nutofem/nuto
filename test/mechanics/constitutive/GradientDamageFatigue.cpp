//
// Created by Thomas Titscher on 10/24/16.
//

#define BOOST_TEST_MODULE GradientDamageFatigueTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "nuto/mechanics/constitutive/laws/GradientDamageFatigueEngineeringStress.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

namespace GradientDamageFatigueTest
{

using namespace NuTo::Constitutive;
BOOST_AUTO_TEST_CASE(AccessParameters)
{
    NuTo::GradientDamageFatigueEngineeringStress law;

    NuTo::ConstitutiveBase& lawBase = law;

    // some parameters from GradientDamage
    lawBase.SetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    lawBase.SetParameterDouble(eConstitutiveParameter::POISSONS_RATIO, .1337);
    lawBase.SetParameterDouble(eConstitutiveParameter::COMPRESSIVE_STRENGTH, 12.);

    BOOST_CHECK_EQUAL(lawBase.GetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS), 30000);
    BOOST_CHECK_EQUAL(lawBase.GetParameterDouble(eConstitutiveParameter::POISSONS_RATIO), .1337);
    BOOST_CHECK_EQUAL(lawBase.GetParameterDouble(eConstitutiveParameter::COMPRESSIVE_STRENGTH), 12.);

    // parameters from GradientDamageFatigue
    // some parameters from GradientDamage
    lawBase.SetParameterDouble(eConstitutiveParameter::ENDURANCE_STRESS, 42.);
    lawBase.SetParameterDouble(eConstitutiveParameter::FATIGUE_PARAMETER, .1337);

    BOOST_CHECK_EQUAL(lawBase.GetParameterDouble(eConstitutiveParameter::ENDURANCE_STRESS), 42.);
    BOOST_CHECK_EQUAL(lawBase.GetParameterDouble(eConstitutiveParameter::FATIGUE_PARAMETER), .1337);

}

BOOST_AUTO_TEST_CASE(UseTheRightKappa)
{
    NuTo::GradientDamageFatigueEngineeringStress law;

}

} // namespace GradientDamageFatigueTest