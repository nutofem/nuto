#include "BoostUnitTest.h"

#include "mechanics/elements/ContinuumElement.h"
#include "mechanics/elements/ElementEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/integrationtypes/IntegrationType1D2NGauss2Ip.h"
#include "mechanics/sections/SectionTruss.h"
#include "mechanics/nodes/NodeDof.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/constitutive/laws/HeatConduction.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "mechanics/elements/ElementOutputBlockMatrixDouble.h"
#include "mechanics/elements/ElementOutputBlockVectorDouble.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"

using namespace NuTo;
BOOST_AUTO_TEST_CASE(check_heat_conduction1D)
{
    Structure structure(1);
    structure.SetNumTimeDerivatives(1);
    std::vector<NodeBase*> nodes;
    std::map<Node::eDof, NuTo::NodeDofInfo> dofInfos;
    dofInfos[Node::eDof::COORDINATES].mDimension = 1;
    dofInfos[Node::eDof::COORDINATES].mNumTimeDerivatives = 1;
    dofInfos[Node::eDof::COORDINATES].mIsDof = false;

    dofInfos[Node::eDof::TEMPERATURE].mDimension = 1;
    dofInfos[Node::eDof::TEMPERATURE].mNumTimeDerivatives = 1;
    dofInfos[Node::eDof::TEMPERATURE].mIsDof = true;

    auto node1 = NodeDof(dofInfos);
    auto node2 = NodeDof(dofInfos);

    Eigen::Matrix<double, 1, 1> coordinate;
    Eigen::Matrix<double, 1, 1> dTemperatureDTime;
    coordinate << 0.0;
    node1.Set(Node::eDof::COORDINATES, 0, coordinate);
    dTemperatureDTime << 1.0;
    node1.Set(Node::eDof::TEMPERATURE, 1, dTemperatureDTime); // set first time derivative for internal gradient calculation

    coordinate << 1.0;
    node2.Set(Node::eDof::COORDINATES, 0, coordinate);
    node2.Set(Node::eDof::TEMPERATURE, 1, dTemperatureDTime);

    nodes.push_back(&node1);
    nodes.push_back(&node2);

    auto truss = Interpolation::eShapeType::TRUSS1D;
    InterpolationType interpolationType(truss, 1);
    interpolationType.AddDofInterpolation(Node::eDof::COORDINATES, Interpolation::eTypeOrder::EQUIDISTANT1);
    interpolationType.AddDofInterpolation(Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);
    IntegrationType1D2NGauss2Ip integrationType;
    interpolationType.UpdateIntegrationType(integrationType);

    ContinuumElement<1> element = ContinuumElement<1>(&structure, nodes, interpolationType);

    ConstitutiveInputMap inputMap;
    std::map<Element::eOutput, std::shared_ptr<ElementOutputBase>> outputMap;
    DofStatus dofStatus;
    std::set<Node::eDof> dofs;
    dofs.insert(Node::eDof::TEMPERATURE);
    dofStatus.SetDofTypes(dofs);
    dofStatus.SetActiveDofTypes(dofs);
    outputMap[Element::eOutput::HESSIAN_0_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(dofStatus);
    outputMap[Element::eOutput::HESSIAN_1_TIME_DERIVATIVE] = std::make_shared<ElementOutputBlockMatrixDouble>(dofStatus);
    outputMap[Element::eOutput::INTERNAL_GRADIENT] = std::make_shared<ElementOutputBlockVectorDouble>(dofStatus);

    SectionTruss area;
    area.SetArea(1.0);
    element.SetSection(area);

    HeatConduction law;
    law.SetParameterDouble(Constitutive::eConstitutiveParameter::THERMAL_CONDUCTIVITY, 1.0);
    law.SetParameterDouble(Constitutive::eConstitutiveParameter::DENSITY, 1.0);
    law.SetParameterDouble(Constitutive::eConstitutiveParameter::HEAT_CAPACITY, 1.0);
    element.SetConstitutiveLaw(law);

    auto interp2 = structure.InterpolationTypeCreate("Truss1D");
    structure.InterpolationTypeAdd(interp2, Node::eDof::TEMPERATURE, Interpolation::eTypeOrder::EQUIDISTANT1);

    element.Evaluate(inputMap, outputMap);

    auto blockhessian0 = outputMap.at(Element::eOutput::HESSIAN_0_TIME_DERIVATIVE)->GetBlockFullMatrixDouble();
    auto hessian0 = blockhessian0.Get("Temperature", "Temperature");
    Eigen::Matrix<double, 2, 2> expected_hessian0;
    expected_hessian0 << 1.0, -1.0, -1.0, 1.0;
    BOOST_CHECK_SMALL((hessian0 - expected_hessian0).norm(), 1e-15);

    auto blockhessian1 = outputMap.at(Element::eOutput::HESSIAN_1_TIME_DERIVATIVE)->GetBlockFullMatrixDouble();
    auto hessian1 = blockhessian1.Get("Temperature", "Temperature");
    Eigen::Matrix<double, 2, 2> expected_hessian1;
    expected_hessian1 << 1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0;
    BOOST_CHECK_SMALL((hessian1 - expected_hessian1).norm(), 1e-15);

    auto blockgradient = outputMap.at(Element::eOutput::INTERNAL_GRADIENT)->GetBlockFullVectorDouble();
    auto gradient = blockgradient.Get("Temperature");
    Eigen::Matrix<double, 2, 1> expected_gradient;
    expected_gradient << 0.5, 0.5;
    BOOST_CHECK_SMALL((gradient - expected_gradient).norm(), 1e-15);
}
