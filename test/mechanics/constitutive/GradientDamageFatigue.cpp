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
#include "../../tools/TypeTraits.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

namespace GradientDamageFatigueTest
{

BOOST_AUTO_TEST_CASE(CopyMove)
{
    NuTo::Test::Copy<NuTo::GradientDamageFatigueEngineeringStress>();
    NuTo::Test::Move<NuTo::GradientDamageFatigueEngineeringStress>();
}


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


NuTo::GradientDamageFatigueEngineeringStress GetLaw()
{
    NuTo::GradientDamageFatigueEngineeringStress law;
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.0);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::ENDURANCE_STRESS, 0.0);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::FATIGUE_PARAMETER, 1.);
    law.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW,
                           static_cast<double>(NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING));
    return law;
}

double GetK0(const NuTo::ConstitutiveBase& rLaw)
{
    return rLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH)
         / rLaw.GetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS);

}

//! @brief this tests checks the case kappa = nonlocaleqstrain
BOOST_AUTO_TEST_CASE(StaticLoading)
{
    auto law = GetLaw();
    auto iplaw = law.CreateIPLaw();
    double k0 = GetK0(iplaw->GetConstitutiveLaw());


    

}

} // namespace GradientDamageFatigueTest