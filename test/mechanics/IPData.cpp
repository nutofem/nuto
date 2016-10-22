//
// Created by Thomas Titscher on 10/22/16.
//

#define BOOST_TEST_MODULE IPDataTest
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include "nuto/mechanics/elements/IPData.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/constitutive/staticData/IPConstitutiveLaw.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/laws/LinearElasticEngineeringStress.h"

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}


const NuTo::LinearElasticEngineeringStress* AsLinearElastic(const NuTo::Constitutive::IPConstitutiveLawBase& rLaw)
{
    const auto* law = dynamic_cast<const NuTo::LinearElasticEngineeringStress*>(&rLaw.GetConstitutiveLaw());
    BOOST_CHECK(law != nullptr);
    return law;
}


BOOST_AUTO_TEST_CASE(IPData_Setup_Test)
{
    NuTo::LinearElasticEngineeringStress law;
    NuTo::IntegrationType1D2NGauss1Ip integrationType1;
    NuTo::IntegrationType1D2NGauss2Ip integrationType2;

    NuTo::IPData data(integrationType1);
    BOOST_CHECK_THROW(data.GetIPConstitutiveLaw(0), NuTo::MechanicsException);
    BOOST_CHECK(not data.HasConstitutiveLawAssigned(0));

    data.SetConstitutiveLaw(law);
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data.GetIPConstitutiveLaw(0)));
    BOOST_CHECK(not data.HasConstitutiveLawAssigned(0));
    BOOST_CHECK(not data.HasConstitutiveLawAssigned(1));

    data.SetIntegrationType(integrationType2);
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data.GetIPConstitutiveLaw(0)));
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data.GetIPConstitutiveLaw(1)));
    BOOST_CHECK(not data.HasConstitutiveLawAssigned(2));

    data.SetIntegrationType(integrationType1);
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data.GetIPConstitutiveLaw(0)));
    BOOST_CHECK(not data.HasConstitutiveLawAssigned(1));
}

BOOST_AUTO_TEST_CASE(IPData_Copy_Move)
{
    NuTo::LinearElasticEngineeringStress law;
    NuTo::IntegrationType1D2NGauss1Ip integrationType;

    NuTo::IPData data(integrationType);
    data.SetConstitutiveLaw(law);

    NuTo::IPData data2(data);   // copy construction
    BOOST_CHECK_EQUAL(&integrationType, &data2.GetIntegrationType());
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data2.GetIPConstitutiveLaw(0)));

    NuTo::IPData data3(integrationType);
    data3 = data;               // copy assignment
    BOOST_CHECK_EQUAL(&integrationType, &data3.GetIntegrationType());
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data3.GetIPConstitutiveLaw(0)));

    NuTo::IPData data4(std::move(data));  // move construction
    BOOST_CHECK_EQUAL(&integrationType, &data4.GetIntegrationType());
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data4.GetIPConstitutiveLaw(0)));

    NuTo::IPData data5(integrationType);
    data5 = std::move(data2);    // move assignment
    BOOST_CHECK_EQUAL(&integrationType, &data5.GetIntegrationType());
    BOOST_CHECK_EQUAL(&law, AsLinearElastic(data5.GetIPConstitutiveLaw(0)));
}
