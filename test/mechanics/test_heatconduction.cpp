#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#define BOOST_TEST_MODULE HeatConductionTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}

BOOST_AUTO_TEST_CASE(tangent_matrix)
{
    // set up constitutive law
    NuTo::HeatConduction heat_conduction;
    std::function<void(double)> setConductivity = std::bind(
            &NuTo::ConstitutiveBase::SetParameterDouble,
            &heat_conduction,
            NuTo::Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 
            std::placeholders::_1);

    double conductivity = 5.0;
    setConductivity(conductivity);

    NuTo::ConstitutiveInputMap input_map;
    NuTo::ConstitutiveOutputMap output_map;
    NuTo::ConstitutiveMatrix<3, 3> tangent;
    output_map[NuTo::Constitutive::Output::D_HEAT_FLUX_D_TEMPERATURE_GRADIENT]
            = &tangent;

    // evaluate tangent
    NuTo::ElementBase* some_element = nullptr;
    heat_conduction.Evaluate<3>(some_element, 0, input_map, output_map);

    // compare to expected output
    Eigen::Matrix<double, 3, 3> matrix;
    matrix << conductivity, 0, 0,
              0, conductivity, 0,
              0, 0, conductivity;
    BOOST_CHECK_EQUAL(tangent, matrix);

    // should throw exception, because conductivity is below zero
    conductivity = -1.0;
    BOOST_CHECK_THROW(setConductivity(conductivity), NuTo::MechanicsException);
}

BOOST_AUTO_TEST_CASE(dofCombinations)
{
    // combinations HeatConduction can compute
    auto rowDof = NuTo::Node::TEMPERATURE;
    auto colDof = NuTo::Node::TEMPERATURE;
    int timeDerivative = 0;
    NuTo::HeatConduction heatConduction;
    BOOST_CHECK(heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));
    timeDerivative = 1;
    BOOST_CHECK(heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));

    // combinations it can't
    colDof = NuTo::Node::DISPLACEMENTS;
    BOOST_CHECK(!heatConduction.CheckDofCombinationComputable(rowDof, colDof, timeDerivative));
}
