#include <stdlib.h>
#include <math.h>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/MechanicsEnums.h"
#include "nuto/mechanics/elements/ElementBase.h"
#include "nuto/visualize/VisualizeEnum.h"

#include "nuto/mechanics/mesh/MeshGenerator.h"

#include "nuto/mechanics/constitutive/laws/AdditiveOutput.h"
#include "nuto/mechanics/constitutive/laws/LinearElasticAnisotropic.h"
#include "nuto/mechanics/constitutive/laws/LinearPiezoelectric.h"
#include "nuto/mechanics/constitutive/laws/LinearDielectric.h"
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOBase.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveScalar.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStrain.h"
#include "nuto/mechanics/constitutive/inputoutput/EngineeringStress.h"
#include "nuto/mechanics/constraints/ConstraintCompanion.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"

#include <boost/filesystem.hpp>
#include "BoostUnitTest.h"


/*
*   Test of a piezoelectric material block
*   uniaxially loaded, bottom electrically grounded
*/

BOOST_AUTO_TEST_CASE(displacementBoundary)
{
    double lX = 0.01;
    double lY = 0.01;
    double lZ = 0.01;
    std::vector<double> blockSize = {lX, lY, lZ};
    std::vector<int> numElements = {2, 2, 2};

    NuTo::Structure s(3);
    std::pair<int, int> gridIds = NuTo::MeshGenerator::Grid(s, blockSize, numElements);

    s.SetNumTimeDerivatives(0);
    s.SetShowTime(false);
    s.SetVerboseLevel(0);

    //    int grp_AllElements = gridIds.first;
    int gridInterpolationType = gridIds.second;
    s.InterpolationTypeAdd(gridInterpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.InterpolationTypeAdd(gridInterpolationType, NuTo::Node::eDof::ELECTRICPOTENTIAL,
                           NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    s.ElementTotalConvertToInterpolationType();

    // ***************************************************************************************************
    //          material data (the famous piezoelectric PZT, strain-charge formulation)
    // ***************************************************************************************************

    double densityPZT = 7600; //[kg/m^2]

    Eigen::MatrixXd compliancePZT(6, 6);
    compliancePZT << 11.5, -3.7, -4.8, 0, 0, 0, -3.7, 11.5, -4.8, 0, 0, 0, -4.8, -4.8, 13.5, 0, 0, 0, 0, 0, 0, 31.9, 0,
            0, 0, 0, 0, 0, 31.9, 0, 0, 0, 0, 0, 0, 30.4;
    compliancePZT *= 1.0e-12; // [m^2 / N]

    Eigen::MatrixXd piezoD_PZT(3, 6);
    piezoD_PZT << 0, 0, 0, 0, 330, 0, 0, 0, 0, 330, 0, 0, -97, -97, 225, 0, 0, 0;
    piezoD_PZT *= 1.0e-12; // [C / N]

    Eigen::MatrixXd dielectric_T_PZT(3, 3);
    dielectric_T_PZT << 1290, 0, 0, 0, 1290, 0, 0, 0, 1000;
    dielectric_T_PZT *= 8.854e-12; // [C/(Vm)]

    // **************************************************************************
    //          material data conversion to stress-charge formulation
    // **************************************************************************

    Eigen::MatrixXd myStiffness(6, 6);
    myStiffness = compliancePZT.inverse();

    Eigen::MatrixXd myPiezo(3, 6);
    myPiezo = piezoD_PZT * myStiffness;

    Eigen::MatrixXd myDielectric(3, 3);
    myDielectric = dielectric_T_PZT - myPiezo * (piezoD_PZT.transpose());

    // **************************************************************************
    //          material law
    // **************************************************************************

    int myLaw = s.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_PIEZOELECTRIC);
    s.ConstitutiveLawSetParameterDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::DENSITY, densityPZT);
    s.ConstitutiveLawSetParameterMatrixDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::STIFFNESS,
                                              myStiffness);
    s.ConstitutiveLawSetParameterMatrixDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::DIELECTRIC_TENSOR,
                                              myDielectric);
    s.ConstitutiveLawSetParameterMatrixDouble(myLaw, NuTo::Constitutive::eConstitutiveParameter::PIEZOELECTRIC_TENSOR,
                                              myPiezo);

    s.ElementTotalSetConstitutiveLaw(myLaw);

    s.NodeBuildGlobalDofs();

    // ***********************************************************
    //          create groups for boundary conditions etc.
    // ***********************************************************

    // bottom
    int grpNodes_Bottom = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Bottom, 2, 0 - 1e-5, 0 + 1e-5);

    // top
    int grpNodes_Top = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Top, 2, lZ - 1e-5, lZ + 1e-5);

    // left
    int grpNodes_Left = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Left, 0, 0 - 1e-5, 0 + 1e-5);

    // right
    int grpNodes_Right = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Right, 0, lX - 1e-5, lX + 1e-5);

    // front
    int grpNodes_Front = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Front, 1, 0 - 1e-5, 0 + 1e-5);

    // back
    int grpNodes_Back = s.GroupCreate(NuTo::eGroupId::Nodes);
    s.GroupAddNodeCoordinateRange(grpNodes_Back, 1, lY - 1e-5, lY + 1e-5);

    // unite all boundary surfaces
    int grpUnionHelper = s.GroupUnion(grpNodes_Top, grpNodes_Bottom);
    grpUnionHelper = s.GroupUnion(grpUnionHelper, grpNodes_Left);
    grpUnionHelper = s.GroupUnion(grpUnionHelper, grpNodes_Right);
    grpUnionHelper = s.GroupUnion(grpUnionHelper, grpNodes_Front);
    int grpNodesBoundary = s.GroupUnion(grpUnionHelper, grpNodes_Back);


    // *************************************************************
    //          create homogeneous deformation test data
    // *************************************************************

    // random values (strain in Voigt notation)
    Eigen::VectorXd testStrain(6);
    testStrain << 1.2, 3.2, 5.4, 1.24, 0.51, 0.86;

    // (cartesian tensor components of strain)
    Eigen::MatrixXd testStrainMatrix(3, 3);
    testStrainMatrix << testStrain(0), 0.5 * testStrain(5), 0.5 * testStrain(4), 0.5 * testStrain(5), testStrain(1),
            0.5 * testStrain(3), 0.5 * testStrain(4), 0.5 * testStrain(3), testStrain(2);

    // random values (electric field)
    Eigen::Vector3d testEField;
    testEField << 0.3545, 0.89514, 1.021564;

    // calculate resulting stresses and electrical displacements (=DField)
    Eigen::VectorXd testStress(6);
    testStress = myStiffness * testStrain - myPiezo.transpose() * testEField;
    Eigen::MatrixXd testStressMatrix(3, 3);
    testStressMatrix << testStress(0), testStress(5), testStress(4), testStress(5), testStress(1), testStress(3),
            testStress(4), testStress(3), testStress(2);

    Eigen::VectorXd testDField(3);
    testDField = myPiezo * testStrain + myDielectric * testEField;

    // *************************************
    //          boundary conditions
    // *************************************

    Eigen::Vector3d directionX;
    directionX << 1., 0., 0.;

    Eigen::Vector3d directionY;
    directionY << 0., 1., 0.0;

    Eigen::Vector3d directionZ;
    directionZ << 0., 0.0, 1.0;

    // electrical conditions at the boundaries are derived from the test data
    std::vector<int> boundaryNodeIds = s.GroupGetMemberIds(grpNodesBoundary);
    for (int nodeId : boundaryNodeIds)
    {
        const auto& curNode = *s.NodeGetNodePtr(nodeId);
        s.Constraints().Add(
                NuTo::Node::eDof::ELECTRICPOTENTIAL,
                NuTo::Constraint::Value(curNode, testEField.dot(curNode.Get(NuTo::Node::eDof::COORDINATES))));
    }

    // mechanical part: either displ. or forces or mixed
    // mechanical conditions at the boundaries are derived from the test data
    for (int nodeId : boundaryNodeIds)
    {
        const auto& curNode = *s.NodeGetNodePtr(nodeId);
        Eigen::VectorXd displ = testStrainMatrix * curNode.Get(NuTo::Node::eDof::COORDINATES);

        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                            NuTo::Constraint::Component(curNode, {NuTo::eDirection::X}, displ[0]));
        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                            NuTo::Constraint::Component(curNode, {NuTo::eDirection::Y}, displ[1]));
        s.Constraints().Add(NuTo::Node::eDof::DISPLACEMENTS,
                            NuTo::Constraint::Component(curNode, {NuTo::eDirection::Z}, displ[2]));
    }

    //    part or all displacement boundary conditions could be replaced with traction boundary conditions below

    //    Eigen::VectorXd loadVector(3);
    //    loadVector << 0, 0, -1;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Bottom,testStressMatrix*loadVector);
    //    loadVector << 0, 0, 1;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Top,testStressMatrix*loadVector);
    //    loadVector << -1, 0, 0;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Left,testStressMatrix*loadVector);
    //    loadVector << 1, 0, 0;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Right,testStressMatrix*loadVector);
    //    loadVector << 0, -1, 0;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Front,testStressMatrix*loadVector);
    //    loadVector << 0, 1, 0;
    //    s.LoadSurfaceConstDirectionCreate3D(0,grp_AllElements,grpNodes_Back,testStressMatrix*loadVector);

    // *************************************
    //          visualize
    // *************************************

    int visualizationGroup = s.GroupCreate(NuTo::eGroupId::Elements);
    s.GroupAddElementsTotal(visualizationGroup);

    //    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ELECTRIC_FIELD);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ELECTRIC_DISPLACEMENT);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ELECTRIC_POTENTIAL);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::DISPLACEMENTS);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    s.AddVisualizationComponent(visualizationGroup, NuTo::eVisualizeWhat::ENGINEERING_STRESS);

    // *************************************
    //          solve static case
    // *************************************

    s.NodeBuildGlobalDofs();

    s.SolveGlobalSystemStaticElastic();

    // *************************************
    //          compare with test case
    // *************************************

    // the result should be constant in all elements (not checked)
    // and also the same for all integration points within one element (also not checked)
    // so just compare data of first element at first integration point to excpected result
    std::vector<int> ElementIds;
    s.ElementGroupGetMembers(visualizationGroup, ElementIds);
    Eigen::MatrixXd resultD =
            s.ElementGetStaticIPData(ElementIds[0], NuTo::IpData::eIpStaticDataType::ELECTRIC_DISPLACEMENT);
    Eigen::MatrixXd resultS =
            s.ElementGetStaticIPData(ElementIds[0], NuTo::IpData::eIpStaticDataType::ENGINEERING_STRAIN);

    BoostUnitTest::CheckEigenMatrix(testDField, resultD.col(0), 1e-6);
    BoostUnitTest::CheckEigenMatrix(testStrain, resultS.col(0), 1e-6);
}
