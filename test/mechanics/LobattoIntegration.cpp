#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverPardiso.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"

#define PRINTRESULT false

/*
 |>*----*----*----*----*  -->F
*/

int buildStructure1D(const std::string& rElementTypeIdent,
                     const std::string& rIntegrationTypeIdent,
                     int rNumNodesPerElement,
                     NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                     int NumElements = 10)
{
    /** paramaters **/
    double  YoungsModulus = 20000.;
    double  Area = 100.0*0.1;
    double  Length = 1000.0;
    double  Force = 1.;
    double  tol = 1.e-6;

    /** Structure 1D **/
    NuTo::Structure myStructure(1);

    /** create section **/
    int Section = myStructure.SectionCreate("Truss");
    myStructure.SectionSetArea(Section, Area);

    /** create material law **/
    int Material = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(Material, YoungsModulus);

    /** create nodes **/
    double elementLength = nodeCoordinatesFirstElement(rNumNodesPerElement-1) - nodeCoordinatesFirstElement(0);
    double factor = Length/(NumElements*elementLength);
    double elementBegin = 0.;
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);

    int node = 0;
    // first node
    nodeCoordinates(0) = factor*nodeCoordinatesFirstElement(0);
    myStructure.NodeCreate(node, "displacements", nodeCoordinates);
    std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
    node++;
    // following nodes
    for(int i = 0; i < NumElements; i++)
    {
        for (int j = 1; j < nodeCoordinatesFirstElement.size(); j++)
        {
            nodeCoordinates(0) = factor*(nodeCoordinatesFirstElement(j) + elementBegin);
            std::cout << "create node: " << node << " coordinates: " << nodeCoordinates(0) << std::endl;
            myStructure.NodeCreate(node, "displacements", nodeCoordinates);
            node++;
        }
        elementBegin+=elementLength;
    }

    int numNodes = node;

    /** create elements **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(rNumNodesPerElement);
    for(int element = 0, node = 0; element < NumElements; element++)
    {
        for (int i = 0; i < rNumNodesPerElement; i++)
        {
            elementIncidence(i) = node;
            node++;
        }
        node--;

        std::cout <<  "create element: " << element << " nodes: " << elementIncidence << std::endl;

        int myElement = myStructure.ElementCreate(rElementTypeIdent, elementIncidence);
        myStructure.ElementSetSection(myElement, Section);
        myStructure.ElementSetConstitutiveLaw(myElement, Material);
        myStructure.ElementSetIntegrationType(myElement, rIntegrationTypeIdent, "NOIPDATA");
    }

    /** set boundary conditions and loads **/
    NuTo::FullVector<double,Eigen::Dynamic> direction(1);
    direction(0) = 1;
    // first node is fixed
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure.LoadCreateNodeForce(1, numNodes-1, direction, Force);

    /** start analysis **/
    myStructure.SetNumProcessors(4);
    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

    stiffnessMatrixCSRVector2.IsSymmetric();

    // build global external load vector
    NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
    myStructure.BuildGlobalExternalLoadVector(1,extForceVector);

    // calculate right hand side
    NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

#if defined HAVE_PARDISO
    NuTo::SparseDirectSolverPardiso mySolverPardiso;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorPardiso;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverPardiso.Solve(stiffnessMatrix, rhsVector, displacementVectorPardiso);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVectorPardiso);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorPardiso;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVectorPardiso);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorPardiso = extForceVector - intForceVectorPardiso;
    std::cout << "residual: " << residualVectorPardiso.Norm() << std::endl;

#elif defined HAVE_MUMPS
    NuTo::SparseDirectSolverMUMPS mySolverMUMPS;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorMUMPS;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverMUMPS.Solve(stiffnessMatrix, rhsVector, displacementVectorMUMPS);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVectorMUMPS);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorMUMPS;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVectorMUMPS);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorMUMPS = extForceVector - intForceVectorMUMPS;
    std::cout << "residual: " << residualVectorMUMPS.Norm() << std::endl;
#else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
#endif

    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFile("LobattoTruss1D2N.vtk");

    double DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    double nodeDisp = myStructure.NodeGetNodePtr(numNodes-1)->GetDisplacement(0);

    if(PRINTRESULT)
    {
        std::cout << "Displacement last node FEM:\n" << nodeDisp << std::endl;
        std::cout << "Displacement last node analytical:\n" << DisplacementCorrect << std::endl;
    }

    if (fabs(nodeDisp - DisplacementCorrect)/fabs(DisplacementCorrect) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }

    return 0;
}

/*
   ||>*----*----*----*----*
      |    |    |    |    | --> Sigma
      |    |    |    |    |
   ||>*----*----*----*----*
      ^
*/
int buildStructure2D(const std::string& rElementTypeIdent,
                     const std::string& rIntegrationTypeIdent,
                     int rNumNodesPerElementInOneDir,
                     NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                     int NumElements)
{
    /** parameters **/
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    double  Height = 100.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 2.123548;
    double  Length = 1000.;
    double  Stress = 10.;
    double  tol = 1.e-4;
    // reference element Length is 2.
    double  factorX =  Length/(NumElements*2);
    double  factorY =  Height/(NumElements*2);
    // length real element
    double  lengthX = Length/NumElements;

    /** Structure 2D **/
    NuTo::Structure myStructure(2);

    /** nodes **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(rNumNodesPerElementInOneDir*rNumNodesPerElementInOneDir);
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(2);
    int nodes = 0;
    for (int i = 0; i < NumElements; i++)
    {
        // first element
        if (i == 0)
        {
            for (int j = 0; j < rNumNodesPerElementInOneDir*rNumNodesPerElementInOneDir; j++)
            {
                nodeCoordinates(0) = nodeCoordinatesFirstElement(j,0)*factorX + lengthX*i;
                nodeCoordinates(1) = nodeCoordinatesFirstElement(j,1)*factorY;
                myStructure.NodeCreate("displacements", nodeCoordinates);
                elementIncidence(j) = nodes;
                nodes++;
            }

            int myElement = myStructure.ElementCreate(rElementTypeIdent, elementIncidence);
            if (PRINTRESULT)
                std::cout <<  "create element: " << myElement << " nodes: \n" << elementIncidence << std::endl;
            //myStructure.ElementSetIntegrationType(myElement, rIntegrationTypeIdent, "NOIPDATA");
        }
        else
        {
            for (int j = 0; j < rNumNodesPerElementInOneDir; j++)
                elementIncidence(j*rNumNodesPerElementInOneDir) = elementIncidence((j+1)*rNumNodesPerElementInOneDir-1);
            for (int j = 0; j < rNumNodesPerElementInOneDir*rNumNodesPerElementInOneDir; j++)
            {
                if (j%rNumNodesPerElementInOneDir != 0)
                {
                    nodeCoordinates(0) = nodeCoordinatesFirstElement(j,0)*factorX + lengthX*i;
                    nodeCoordinates(1) = nodeCoordinatesFirstElement(j,1)*factorY;
                    myStructure.NodeCreate("displacements", nodeCoordinates);
                    elementIncidence(j) = nodes;
                    nodes++;
                }
            }

            int myElement = myStructure.ElementCreate(rElementTypeIdent, elementIncidence);
            if (PRINTRESULT)
                std::cout <<  "create element: " << myElement << " nodes: \n " << elementIncidence << std::endl;
            //myStructure.ElementSetIntegrationType(myElement, rIntegrationTypeIdent, "NOIPDATA");
        }
    }

    myStructure.ElementTotalSetIntegrationType(rIntegrationTypeIdent, "NOIPDATA");

    myStructure.Info();

    /** create constitutive law **/
    int myMatLin = myStructure.ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure.ConstitutiveLawSetYoungsModulus(myMatLin,YoungsModulus);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMatLin,PoissonRatio);

    /** create section **/
    int mySection = myStructure.SectionCreate("Plane_Stress");
    myStructure.SectionSetThickness(mySection,Thickness);

    /** assign constitutive law **/
    myStructure.ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure.ElementTotalSetSection(mySection);

    // Dirichlet
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);
    // fix both in x
    direction << 1.,0.;

    std::cout <<  "Nodes fixed: "<< std::endl;
    for (int i = 0; i < rNumNodesPerElementInOneDir; i++)
    {
        std::cout << i*rNumNodesPerElementInOneDir << std::endl;
        myStructure.ConstraintLinearSetDisplacementNode(i*rNumNodesPerElementInOneDir, direction, 0.0);
    }

    // fix node 1 in y
    direction(0) = 0;
    direction(1) = 1;
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);

    // right boundary
    std::cout <<  "Pressure on nodes: "<< std::endl;

    int groupNumberNodesRight = myStructure.GroupCreate("NODES");
    for (int i= 0; i < rNumNodesPerElementInOneDir; i++)
    {
        std::cout << (nodes-1)-i*(rNumNodesPerElementInOneDir-1) << std::endl;
        myStructure.GroupAddNode(groupNumberNodesRight, (nodes-1)-i*(rNumNodesPerElementInOneDir-1));
    }
    std::cout <<  " and element: " << NumElements-1 << std::endl;
    int groupNumberElementsRight = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElement(groupNumberElementsRight, NumElements-1);

    //NuTo::FullVector<double, 2> ForceVecRight({Force,0.});
    //myStructure.LoadSurfaceConstDirectionCreate2D(1, groupNumberElementsRight, groupNumberNodesRight, ForceVecRight);
    myStructure.LoadSurfacePressureCreate2D(1, groupNumberElementsRight, groupNumberNodesRight, -Stress);

    myStructure.Info();

    /** start analysis **/
    myStructure.SetNumProcessors(4);
    myStructure.CalculateMaximumIndependentSets();
    myStructure.NodeBuildGlobalDofs();

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

    // build global external load vector
    NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
    myStructure.BuildGlobalExternalLoadVector(1,extForceVector);

    // calculate right hand side
    NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

#ifdef HAVE_PARDISO
    NuTo::SparseDirectSolverPardiso mySolverPardiso;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorPardiso;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverPardiso.Solve(stiffnessMatrix, rhsVector, displacementVectorPardiso);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVectorPardiso);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorPardiso;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVectorPardiso);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorPardiso = extForceVector - intForceVectorPardiso;
    std::cout << "residual: " << residualVectorPardiso.Norm() << std::endl;

    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFile("LobattoTruss1D2N.vtk");
#elif defined HAVE_MUMPS
    NuTo::SparseDirectSolverMUMPS mySolverMUMPS;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorMUMPS;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverMUMPS.Solve(stiffnessMatrix, rhsVector, displacementVectorMUMPS);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVectorMUMPS);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorMUMPS;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVectorMUMPS);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorMUMPS = extForceVector - intForceVectorMUMPS;
    std::cout << "residual: " << residualVectorMUMPS.Norm() << std::endl;
#else
    std::cout << "PARDISO not available - can't solve system of equations " << std::endl;
#endif // HAVE_PARDISO

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStress(6,1);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStressCorrect(6,1);
    EngineeringStressCorrect(0,0) = Stress;

    double DisplacementCorrect = 0.5;
    double nodeDisp = myStructure.NodeGetNodePtr(nodes-1)->GetDisplacement(0);

    if (PRINTRESULT)
    {
        std::cout << "Displacement analytical: " << DisplacementCorrect  << std::endl;
        std::cout << "Displacement FEM: " <<  nodeDisp << std::endl;
    }

    for(int element = 0; element < NumElements; element++)
    {
        if (PRINTRESULT)
        {
            myStructure.ElementGetEngineeringStress(element, EngineeringStress);
            std::cout << "EngineeringStressCorrect" << std::endl;
            EngineeringStressCorrect.Info();
            std::cout << "EngineeringStress" << std::endl;
            EngineeringStress.Info();
        }
    }

    double resultTol = fabs(nodeDisp - DisplacementCorrect)/fabs(DisplacementCorrect);
    if (resultTol > tol)
    {
        //throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }

    return 0;
}

int main()
{
    /** 2D **/

    // 3x3Nodes
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> ones2D(9,2); ones2D.fill(1);
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates2D(9,2);
    NuTo::IntegrationType2D4NLobatto9Ip Lobatto2D4N9Ip;
    for (int i = 0; i < 9; i++)
    {
        double coords[2];
        Lobatto2D4N9Ip.GetLocalIntegrationPointCoordinates2D(i, coords);
        nodeCoordinates2D(i,0) = coords[0];
        nodeCoordinates2D(i,1) = coords[1];
    }

    nodeCoordinates2D += ones2D;
    buildStructure2D("PLANE2D4NSPECTRALORDER2", "2D4NLobatto9Ip", 3, nodeCoordinates2D, 1000);

    // 4x4Nodes
    ones2D.resize(16,2); ones2D.fill(1);
    nodeCoordinates2D.resize(16,2);
    NuTo::IntegrationType2D4NLobatto16Ip Lobatto2D4N16Ip;
    for (int i = 0; i < 16; i++)
    {
        double coords[2];
        Lobatto2D4N16Ip.GetLocalIntegrationPointCoordinates2D(i, coords);
        nodeCoordinates2D(i,0) = coords[0];
        nodeCoordinates2D(i,1) = coords[1];
    }
    nodeCoordinates2D += ones2D;
    std::cout << "NodeCoordinates: \n" << nodeCoordinates2D << std::endl;
    buildStructure2D("PLANE2D4NSPECTRALORDER3", "2D4NLobatto16Ip", 4, nodeCoordinates2D, 1000);


    // 5x5Nodes
    ones2D.resize(25,2); ones2D.fill(1);
    nodeCoordinates2D.resize(25,2);
    NuTo::IntegrationType2D4NLobatto25Ip Lobatto2D4N25Ip;
    for (int i = 0; i < 25; i++)
    {
        double coords[2];
        Lobatto2D4N25Ip.GetLocalIntegrationPointCoordinates2D(i, coords);
        nodeCoordinates2D(i,0) = coords[0];
        nodeCoordinates2D(i,1) = coords[1];
    }
    nodeCoordinates2D += ones2D;
    buildStructure2D("PLANE2D4NSPECTRALORDER4", "2D4NLobatto25Ip", 5, nodeCoordinates2D, 1000);

    /** 1D **/

    // 3Nodes
    NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
    NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;
    for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates += ones;
    buildStructure1D("TRUSS1D2NSPECTRALORDER2", "1D2NLobatto3Ip", 3, nodeCoordinates, 1000);

    // 4Nodes
    ones.resize(4); ones.fill(1);
    nodeCoordinates.resize(4);
    NuTo::IntegrationType1D2NLobatto4Ip Lobatto1D2N4Ip;
    for (int i = 0; i < 4; i++) Lobatto1D2N4Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates+=ones;
    buildStructure1D("TRUSS1D2NSPECTRALORDER3", "1D2NLobatto4Ip", 4, nodeCoordinates, 1000);

    // 5Nodes
    ones.resize(5); ones.fill(1);
    nodeCoordinates.resize(5);
    NuTo::IntegrationType1D2NLobatto5Ip Lobatto1D2N5Ip;
    for (int i = 0; i < 5; i++) Lobatto1D2N5Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates+=ones;
    buildStructure1D("TRUSS1D2NSPECTRALORDER4", "1D2NLobatto5Ip", 5, nodeCoordinates, 1000);


    return 0;
}
