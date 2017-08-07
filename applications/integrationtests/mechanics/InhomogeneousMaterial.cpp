#include <iostream>

#include "mechanics/constitutive/laws/LinearElasticInhomogeneous.h"

#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/Assembler.h"

#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/MechanicsEnums.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/groups/Group.h"
#include "mechanics/sections/SectionTruss.h"
#include "BoostUnitTest.h"

/*
 * 1D Unixaxial tension test, prescribed displacement at boundaries
 * Youngs modulus spatially varying
 */
BOOST_AUTO_TEST_CASE(InhomogeneousMaterial1D)
{
    int numElm = 10;
    double lX = 1.0;
    double E = 1.2;
    double rightBoundaryValue = 0.9;

    NuTo::Structure s(1);
    std::pair<int, int> gridIds = NuTo::MeshGenerator::Grid(s, {lX}, {numElm} );

    int grp_AllNodes = s.GroupCreateNodeGroupFromElements(gridIds.first);
    s.InterpolationTypeAdd(gridIds.second,NuTo::Node::eDof::DISPLACEMENTS,NuTo::Interpolation::eTypeOrder::LOBATTO2);
    s.ElementTotalConvertToInterpolationType();

    int myLaw = 0;
    s.ConstitutiveLawCreate(myLaw, NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_INHOMOGENEOUS);
    s.ConstitutiveLawSetParameterDouble(myLaw,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, E);
    s.ConstitutiveLawSetParameterDouble(myLaw,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,0);
    s.ConstitutiveLawSetParameterDouble(myLaw,NuTo::Constitutive::eConstitutiveParameter::DENSITY, 1);

    NuTo::LinearElasticInhomogeneous* myLawPtr =
            static_cast<NuTo::LinearElasticInhomogeneous*>(s.ConstitutiveLawGetConstitutiveLawPtr(myLaw));

    // Varying Youngs modulus
    std::function<double(Eigen::VectorXd)> functionYoungsModulus = [E,lX](Eigen::VectorXd x)
    {
        return E*1/(x(0)/lX+1.);
    };
    myLawPtr->SetYoungsModulus(functionYoungsModulus);

    s.ElementTotalSetConstitutiveLaw(myLaw);

    s.NodeBuildGlobalDofs();

    NuTo::Group<NuTo::NodeBase> grpNodes_Left = s.GroupGetNodeCoordinateRange(NuTo::eDirection::X, 0 - 1e-5, 0 + 1e-5);
    NuTo::Group<NuTo::NodeBase> grpNodes_Right = s.GroupGetNodeCoordinateRange(NuTo::eDirection::X, lX - 1e-5, lX + 1e-5);

    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                NuTo::Constraint::Component(grpNodes_Left, {NuTo::eDirection::X},0.0));

    s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                NuTo::Constraint::Component(grpNodes_Right, {NuTo::eDirection::X},rightBoundaryValue));

    s.ElementTotalSetSection(NuTo::SectionTruss::Create(1.0));


    s.SolveGlobalSystemStaticElastic();

    Eigen::MatrixXd displResult;
    Eigen::MatrixXd coordResult;
    s.NodeGroupGetDisplacements(grp_AllNodes,displResult);
    s.NodeGroupGetCoordinates(grp_AllNodes,coordResult);

    // expected result is C1 * int_0^lX 1/E(x) dx, constant C1 fulfills boundary condition on the right

    Eigen::MatrixXd displExpected = 0. * displResult;
    for (int i=0; i<displResult.size(); i++) {
        double x = coordResult(i) / lX;
        displExpected(i) = rightBoundaryValue * 2 / 3 * (1./2.*x*x + x);
    }

    BoostUnitTest::CheckEigenMatrix(displResult,displExpected);
}
