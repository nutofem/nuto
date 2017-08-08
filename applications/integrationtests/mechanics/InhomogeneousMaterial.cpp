#include <iostream>

#include "mechanics/constitutive/laws/LinearElasticInhomogeneous.h"

#include "mechanics/structures/unstructured/Structure.h"

#include "BoostUnitTest.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/sections/SectionTruss.h"

/*
 * 1D Unixaxial tension test, prescribed displacement at boundaries
 * Youngs modulus spatially varying
 */
BOOST_AUTO_TEST_CASE(InhomogeneousMaterial1D) {
  int numElm = 10;
  double lX = 1.0;
  double E = 1.2;
  double rightBoundaryValue = 0.9;

  NuTo::Structure s(1);
  std::pair<int, int> gridIds = NuTo::MeshGenerator::Grid(s, {lX}, {numElm});

  int grp_AllNodes = s.GroupCreateNodeGroupFromElements(gridIds.first);
  s.InterpolationTypeAdd(gridIds.second, NuTo::Node::eDof::DISPLACEMENTS,
                         NuTo::Interpolation::eTypeOrder::LOBATTO2);
  s.ElementTotalConvertToInterpolationType();

  // Create Law
  int myLawId = s.ConstitutiveLawCreate(
      NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_INHOMOGENEOUS);
  NuTo::LinearElasticInhomogeneous &myLaw =
      *(dynamic_cast<NuTo::LinearElasticInhomogeneous *>(
          s.ConstitutiveLawGetConstitutiveLawPtr(myLawId)));
  myLaw.SetParameterDouble(
      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, 0.);
  myLaw.SetParameterDouble(NuTo::Constitutive::eConstitutiveParameter::DENSITY,
                           1.);
  // Varying Youngs modulus
  std::function<double(Eigen::VectorXd)> functionYoungsModulus =
      [E, lX](Eigen::VectorXd x) { return E * 1 / (x(0) / lX + 1.); };
  myLaw.SetYoungsModulus(functionYoungsModulus);

  s.ElementTotalSetConstitutiveLaw(myLawId);

  s.NodeBuildGlobalDofs();

  NuTo::NodeBase &node_Left = s.NodeGetAtCoordinate(0.);
  NuTo::NodeBase &node_Right = s.NodeGetAtCoordinate(lX);

  s.Constraints().Add(
      NuTo::Node::eDof::DISPLACEMENTS,
      NuTo::Constraint::Component(node_Left, {NuTo::eDirection::X}, 0.0));

  s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                      NuTo::Constraint::Component(node_Right,
                                                  {NuTo::eDirection::X},
                                                  rightBoundaryValue));

  s.ElementTotalSetSection(NuTo::SectionTruss::Create(1.0));

  s.SolveGlobalSystemStaticElastic();

  Eigen::MatrixXd displResult;
  Eigen::MatrixXd coordResult;
  s.NodeGroupGetDisplacements(grp_AllNodes, displResult);
  s.NodeGroupGetCoordinates(grp_AllNodes, coordResult);

  // expected result is C1 * int_0^lX 1/E(x) dx, constant C1 fulfills boundary
  // condition on the right

  Eigen::MatrixXd displExpected = 0. * displResult;
  for (int i = 0; i < displResult.size(); i++) {
    double x = coordResult(i) / lX;
    displExpected(i) = rightBoundaryValue * 2 / 3 * (1. / 2. * x * x + x);
  }

  BoostUnitTest::CheckEigenMatrix(displResult, displExpected);
}
