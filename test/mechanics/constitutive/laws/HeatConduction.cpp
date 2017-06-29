#include "BoostUnitTest.h"

#include "mechanics/constitutive/laws/HeatConduction.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveMatrix.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/nodes/NodeEnum.h"

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
    heat_conduction.Evaluate<3>(input_map, output_map);

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
    BOOST_CHECK_THROW(setConductivity(conductivity), NuTo::Exception);
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
