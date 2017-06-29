#include "BoostUnitTest.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/structures/unstructured/Structure.h"
#include "base/Exception.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constraints/ConstraintCompanion.h"

void CoefficientCheckLinearElasticTriangle(NuTo::Interpolation::eTypeOrder rTypeOrder)
{
    //create structure
    NuTo::Structure myStructure(2);
    //create nodes
    Eigen::MatrixXd nodeCoordinates(2, 8);
    nodeCoordinates <<
            0, 10, 2, 8, 4, 8, 0, 10,
            0, 0, 2, 3, 7, 7, 10, 10;
    myStructure.NodesCreate(nodeCoordinates);

    int interpolationType = myStructure.InterpolationTypeCreate("Triangle2D");
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, rTypeOrder);

    //create element
    Eigen::MatrixXi nodeNumbers(3, 10);
    nodeNumbers <<
            0, 0, 0, 1, 2, 2, 3, 3, 4, 5,
            1, 2, 3, 7, 4, 3, 5, 7, 5, 7,
            3, 6, 2, 3, 6, 4, 4, 5, 6, 6;

    myStructure.ElementsCreate(interpolationType, nodeNumbers);
    myStructure.ElementTotalConvertToInterpolationType();


    //Calculate maximum independent sets for parallelization (openmp)
    myStructure.CalculateMaximumIndependentSets();

    //create constitutive law
    int myMatLin = myStructure.ConstitutiveLawCreate("Linear_Elastic_Engineering_Stress");
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, 10);
    myStructure.ConstitutiveLawSetParameterDouble(myMatLin, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.2);

    //create section
    auto mySection = NuTo::SectionPlane::Create(1.0, true);

    //assign constitutive law
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    // add a rather random constraint to get some dependent dofs
    const auto& node = *myStructure.NodeGetNodePtr(0);
    myStructure.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Constraint::Component(node, {NuTo::eDirection::X}));
    myStructure.CalculateMaximumIndependentSets();

    BOOST_CHECK(myStructure.CheckHessian0(1.e-6, 1.e-4));
    BOOST_CHECK(myStructure.ElementCheckHessian0(1.e-6, 1.e-4));
}

BOOST_AUTO_TEST_CASE(CoeffsOrder1)
{
        CoefficientCheckLinearElasticTriangle(NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
}
BOOST_AUTO_TEST_CASE(CoeffsOrder2)
{
        CoefficientCheckLinearElasticTriangle(NuTo::Interpolation::eTypeOrder::EQUIDISTANT2);
}
