#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"

#define BOOST_TEST_MODULE HeatConductionTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

using namespace NuTo;
BOOST_AUTO_TEST_CASE(tangent_matrix)
{
    // set up constitutive law
    HeatConduction heat_conduction;
    std::function<void(double)> setConductivity = std::bind(
            &NuTo::ConstitutiveBase::SetParameterDouble,
            &heat_conduction,
            Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY,
            std::placeholders::_1);

    double conductivity = 5.0;
    setConductivity(conductivity);

    ConstitutiveInputMap input_map;
    ConstitutiveOutputMap output_map;
    output_map[Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT] =
            ConstitutiveIOBase::makeConstitutiveIO<3>(Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT);

    // evaluate tangent
    ElementBase* some_element = nullptr;
    heat_conduction.Evaluate<3>(some_element, 0, input_map, output_map);

    // compare to expected output
    Eigen::Matrix<double, 3, 3> expected_conductivity;
    expected_conductivity << conductivity, 0, 0,
                             0, conductivity, 0,
                             0, 0, conductivity;
    Eigen::Matrix<double, 3, 3> calculated_conductivity =
            *static_cast<ConstitutiveMatrix<3,3>*>(output_map.at(Constitutive::eOutput::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT).get());
    BOOST_CHECK_EQUAL(calculated_conductivity, expected_conductivity);

    // should throw exception, because conductivity is below zero
    conductivity = -1.0;
    BOOST_CHECK_THROW(setConductivity(conductivity), NuTo::MechanicsException);
}

BOOST_AUTO_TEST_CASE(dofCombinations)
{
    // combinations HeatConduction can compute
    auto rowDof = Node::eDof::TEMPERATURE;
    auto colDof = Node::eDof::TEMPERATURE;
    int timeDerivative = 0;
    HeatConduction heatConduction;
    BOOST_CHECK(heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));
    timeDerivative = 1;
    BOOST_CHECK(heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));

    // combinations it can't
    colDof = NuTo::Node::eDof::DISPLACEMENTS;
    BOOST_CHECK(!heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));
}
