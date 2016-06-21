#include "nuto/mechanics/elements/ContinuumElement.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NGauss1Ip.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include "nuto/mechanics/nodes/NodeDof.h"
#include "nuto/mechanics/constitutive/laws/HeatConduction.h"
#include "nuto/mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "nuto/mechanics/dofSubMatrixStorage/DofStatus.h"
#define BOOST_TEST_MODULE ContinuumElementTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

using namespace NuTo;

// necessary to build with clang when boost has been compiled by gcc
std::string boost::unit_test::ut_detail::normalize_test_case_name(const_string name)
{
    return (name[0] == '&' ? std::string(name.begin()+1, name.size()-1) : std::string(name.begin(), name.size() ));
}


BOOST_AUTO_TEST_CASE(check_evaluate)
{
    Structure structure(1);
    std::vector<NodeBase*> nodes;
    NuTo::NodeDof<1, 0, 0, 0, 1, 0, 0, 0, 0, 0> node1;
    NuTo::NodeDof<1, 0, 0, 0, 1, 0, 0, 0, 0, 0> node2;

    Eigen::Matrix<double, 1, 1> coordinate;
    coordinate << 0.0;
    node1.SetCoordinates1D(coordinate);
    coordinate << 1.0;
    node2.SetCoordinates1D(coordinate);
    nodes.push_back(&node1);
    nodes.push_back(&node2);

    auto elementDataType = ElementData::eElementDataType::CONSTITUTIVELAWIP;
    auto ipDataType = IpData::eIpDataType::NOIPDATA;

    auto truss = Interpolation::eShapeType::TRUSS1D;
    InterpolationType interpolationType(truss, 1);
    interpolationType.AddDofInterpolation(Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    interpolationType.AddDofInterpolation(Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);
    IntegrationType1D2NGauss1Ip integrationType;
    interpolationType.UpdateIntegrationType(integrationType);

    ContinuumElement<1> element = ContinuumElement<1>(&structure,
            nodes, elementDataType, ipDataType, &interpolationType);

    ConstitutiveInputMap inputMap;
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> outputMap;
    DofStatus dofStatus;
    std::set<Node::eDof> dofs;
    dofs.insert(Node::eDof::TEMPERATURE);
    dofStatus.SetDofTypes(dofs);
    dofStatus.SetActiveDofTypes(dofs);
    outputMap[Element::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(dofStatus);

    SectionTruss area;
    area.SetArea(1.0);
    element.SetSection(&area);

    HeatConduction law;
    law.SetParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.0);
    element.SetConstitutiveLaw(&law);

    auto interp2 = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(interp2, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);

    element.Evaluate(inputMap, outputMap);

    auto blockhessian = outputMap[Element::HESSIAN_0_TIME_DERIVATIVE]->GetBlockFullMatrixDouble();
    auto hessian = blockhessian.Get("Temperature", "Temperature");

    Eigen::Matrix<double, 2, 2> expected_hessian;
    expected_hessian << 1.0, -1.0, -1.0, 1.0;
    BOOST_CHECK_EQUAL(hessian, expected_hessian);
}
