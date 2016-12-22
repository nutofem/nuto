//
// Created by Thomas Titscher on 10/24/16.
//

#define BOOST_TEST_MODULE GradientDamageFatigueTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "mechanics/constitutive/laws/GradientDamageFatigueEngineeringStress.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "../../tools/TypeTraits.h"
#include "ConstitutiveTangentTester.h"

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


using namespace NuTo;
using namespace NuTo::Constitutive;

BOOST_AUTO_TEST_CASE(AccessParameters)
{
    GradientDamageFatigueEngineeringStress law;

    ConstitutiveBase& lawBase = law;

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
    law.SetParameterDouble(eConstitutiveParameter::DENSITY, 1.0);
    law.SetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS, 30000);
    law.SetParameterDouble(eConstitutiveParameter::POISSONS_RATIO, 0.0);
    law.SetParameterDouble(eConstitutiveParameter::NONLOCAL_RADIUS, 1);
    law.SetParameterDouble(eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 0.);
    law.SetParameterDouble(eConstitutiveParameter::TENSILE_STRENGTH, 4.);
    law.SetParameterDouble(eConstitutiveParameter::COMPRESSIVE_STRENGTH, 4. * 10);
    law.SetParameterDouble(eConstitutiveParameter::FRACTURE_ENERGY, 0.021);
    law.SetParameterDouble(eConstitutiveParameter::ENDURANCE_STRESS, 0.0);
    law.SetParameterDouble(eConstitutiveParameter::FATIGUE_PARAMETER, 1.);
    law.SetParameterDouble(eConstitutiveParameter::DAMAGE_LAW,
                           static_cast<double>(eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING));
    return law;
}

double GetK0(const ConstitutiveBase& rLaw)
{
    return rLaw.GetParameterDouble(eConstitutiveParameter::TENSILE_STRENGTH)
         / rLaw.GetParameterDouble(eConstitutiveParameter::YOUNGS_MODULUS);

}

template <int TDim>
ConstitutiveInputMap GetInputMap(EngineeringStrain<TDim> rStrain, double rEeq)
{
    ConstitutiveInputMap input;
    input.Add<TDim>(eInput::NONLOCAL_EQ_STRAIN);
    input.at(eInput::NONLOCAL_EQ_STRAIN)->operator[](0) = rEeq;

    input.Add<TDim>(eInput::ENGINEERING_STRAIN);
    input.at(eInput::ENGINEERING_STRAIN)->AsEngineeringStrain<TDim>() = rStrain;

    input.Add<TDim>(eInput::PLANE_STATE);
    return input;
}

//! @brief this tests checks the case kappa = nonlocaleqstrain
BOOST_AUTO_TEST_CASE(StaticLoading)
{
    auto law = GetLaw();
    IPConstitutiveLaw<GradientDamageFatigueEngineeringStress> iplaw(law, Eigen::Vector2d(0.,0.));
    auto& staticData = iplaw.GetStaticData();
    staticData.SetData(Eigen::Vector2d(0.,0.));

    double k0 = GetK0(law);

    auto input = GetInputMap<2>({0.5*k0, 0, 0}, 0.5*k0);
    NuTo::Test::ConstitutiveTangentTester<2> tester(iplaw, 1.e-8, 1.e-6);

    BOOST_CHECK(tester.CheckTangent(input,
                                    eInput::ENGINEERING_STRAIN,
                                    eOutput::ENGINEERING_STRESS,
                                    eOutput::D_ENGINEERING_STRESS_D_ENGINEERING_STRAIN));

    BOOST_CHECK(tester.CheckTangent(input,
                                    eInput::NONLOCAL_EQ_STRAIN,
                                    eOutput::ENGINEERING_STRESS,
                                    eOutput::D_ENGINEERING_STRESS_D_NONLOCAL_EQ_STRAIN));

    BOOST_CHECK(tester.CheckTangent(input,
                                    eInput::ENGINEERING_STRAIN,
                                    eOutput::LOCAL_EQ_STRAIN,
                                    eOutput::D_LOCAL_EQ_STRAIN_D_STRAIN));

}

} // namespace GradientDamageFatigueTest