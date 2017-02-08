#include "BoostUnitTest.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"

using namespace NuTo;
BOOST_AUTO_TEST_CASE(modify_value_in_map)
{
    using namespace NuTo::Constitutive;
    ConstitutiveInputMap inputMap;
    inputMap[eInput::TEMPERATURE] = ConstitutiveIOBase::makeConstitutiveIO<1>(eInput::TEMPERATURE);
    double& innertemp = (*inputMap.at(eInput::TEMPERATURE))[0];
    innertemp = 5.0;
    auto value_in_map = (*inputMap[eInput::TEMPERATURE])[0];
    BOOST_CHECK_EQUAL(value_in_map, 5.0);
}

BOOST_AUTO_TEST_CASE(create_input)
{
    ConstitutiveInputMap inputMap;
    // Scalar
    auto inputType = NuTo::Constitutive::eInput::TEMPERATURE;
    inputMap[inputType] = ConstitutiveIOBase::makeConstitutiveIO<1>(inputType);
    double& innertemp = (*inputMap.at(inputType))[0];
    innertemp = 5.0;
    auto value_in_map = (*inputMap.at(inputType))[0];
    BOOST_CHECK_EQUAL(value_in_map, 5.0);

    // Vector
    inputType = NuTo::Constitutive::eInput::TEMPERATURE_GRADIENT;
    inputMap[inputType] = ConstitutiveIOBase::makeConstitutiveIO<2>(inputType);
    auto& innergrad = static_cast<ConstitutiveVector<2>*>(inputMap.at(inputType).get())->AsVector();
    innergrad[0] = 5.0; innergrad[1] = 7.0;
    double first = (*inputMap.at(inputType))[0];
    double second = (*inputMap.at(inputType))[1];
    BOOST_CHECK_EQUAL(first, 5.0);
    BOOST_CHECK_EQUAL(second, 7.0);

    // Engineering Strain
    inputType = NuTo::Constitutive::eInput::ENGINEERING_STRAIN;
    inputMap[inputType] = ConstitutiveIOBase::makeConstitutiveIO<3>(inputType);
    auto& innerstrain = static_cast<EngineeringStrain<3>*>(inputMap.at(inputType).get())->AsEngineeringStrain3D();
    innerstrain[0] = 1.0; innerstrain[3] = 6.0; innerstrain[5] = 3.0;
    first = (*inputMap.at(inputType))[0];
    second = (*inputMap.at(inputType))[3];
    double third = (*inputMap.at(inputType))[5];
    BOOST_CHECK_EQUAL(first, 1.0);
    BOOST_CHECK_EQUAL(second, 6.0);
    BOOST_CHECK_EQUAL(third, 3.0);
}

BOOST_AUTO_TEST_CASE(copy_construction)
{
    ConstitutiveInputMap inputMap;
    auto inputType = NuTo::Constitutive::eInput::TEMPERATURE;
    inputMap[inputType] = ConstitutiveIOBase::makeConstitutiveIO<1>(inputType);
    double& temp = (*inputMap.at(inputType))[0];
    temp = 5.0;

    ConstitutiveInputMap emptyMap = inputMap;
    BOOST_CHECK_EQUAL((*inputMap.at(inputType))[0], (*emptyMap.at(inputType))[0]);
}
