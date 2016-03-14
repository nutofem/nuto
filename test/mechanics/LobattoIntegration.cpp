#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include <eigen3/Eigen/Core>

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseMatrixCSRSymmetric.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverPardiso.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto3Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto4Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType1D2NLobatto5Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto9Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto16Ip.h"
#include "nuto/mechanics/integrationtypes/IntegrationType2D4NLobatto25Ip.h"

#include "nuto/mechanics/timeIntegration/RungeKutta4.h"

#define PRINTRESULT false

/*
 |>*----*----*----*----*  -->F
*/

NuTo::Structure* buildStructure1D(NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                                  int rNumNodesPerElement,
                                  NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                                  int NumElements,
                                  int timeDers,
                                  double& DisplacementCorrect)
{
    /** paramaters **/
    double  YoungsModulus = 20000.;
    double  Area = 100.0*0.1;
    double  Length = 1000.0;
    double  Force = 1.;

    DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    /** Structure 1D **/
    NuTo::Structure* myStructure = new NuTo::Structure(1);
    myStructure->SetNumTimeDerivatives(timeDers);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    /** create section **/
    int Section = myStructure->SectionCreate("Truss");
    myStructure->SectionSetArea(Section, Area);

    /** create material law **/
    int Material = myStructure->ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure->ConstitutiveLawSetParameterDouble(Material, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, YoungsModulus);

    /** create nodes **/
    double elementLength = nodeCoordinatesFirstElement(rNumNodesPerElement-1) - nodeCoordinatesFirstElement(0);
    double factor = Length/(NumElements*elementLength);
    double elementBegin = 0.;
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);

    int node = 0;
    // first node
    nodeCoordinates(0) = factor*nodeCoordinatesFirstElement(0);
    myStructure->NodeCreate(node, nodeCoordinates);
    if(PRINTRESULT) std::cout << "create node: " << node << " coordinates: " << nodeCoordinates.at(0,0) << std::endl;
    node++;
    // following nodes
    for(int i = 0; i < NumElements; i++)
    {
        for (int j = 1; j < nodeCoordinatesFirstElement.size(); j++)
        {
            nodeCoordinates(0) = factor*(nodeCoordinatesFirstElement(j) + elementBegin);
            if(PRINTRESULT) std::cout << "create node: " << node << " coordinates: " << nodeCoordinates.at(0,0) << std::endl;
            myStructure->NodeCreate(node, nodeCoordinates);
            node++;
        }
        elementBegin+=elementLength;
    }

    int interpolationType = myStructure->InterpolationTypeCreate("TRUSS1D");
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, rElementTypeIdent);
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, rElementTypeIdent);

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

        if(PRINTRESULT) std::cout <<  "create element: " << element << " nodes: " << elementIncidence << std::endl;

        myStructure->ElementCreate(interpolationType, elementIncidence);
    }

    myStructure->ElementTotalConvertToInterpolationType(1e-6,3);
    myStructure->ElementTotalSetSection(Section);
    myStructure->ElementTotalSetConstitutiveLaw(Material);

    /** set boundary conditions and loads **/
    NuTo::FullVector<double,Eigen::Dynamic> direction(1);
    direction(0) = 1;
    // first node is fixed
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure->SetNumLoadCases(1);
    myStructure->LoadCreateNodeForce(0, numNodes-1, direction, Force);

    /** start analysis **/

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    return myStructure;
}


/*
   ||>*----*----*----*----*
      |    |    |    |    | -->
      |    |    |    |    |
   ||>*----*----*----*----*     Sigma
      |    |    |    |    | -->
      |    |    |    |    |
   ||>*----*----*----*----*
      ^
*/
NuTo::Structure* buildStructure2D(NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                     //NuTo::IntegrationType::eIntegrationType rIntegrationTypeIdent,
                     int rNumNodesPerElementInOneDir,
                     NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                     int NumElementsX, int NumElementsY, double& DisplacementCorrect)
{
    /** parameters **/
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    double  Height = 5.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 2.123548;
    double  Length = 10;
    double  Stress = 10.;
    double elementSize = nodeCoordinatesFirstElement(rNumNodesPerElementInOneDir-1) - nodeCoordinatesFirstElement(0);
    // reference element Length is 2.
    double  factorX =  Length/(NumElementsX*elementSize);
    double  factorY =  Height/(NumElementsY*elementSize);

    DisplacementCorrect = (Stress*Length)/YoungsModulus;

    /** Structure 2D **/
    NuTo::Structure* myStructure = new NuTo::Structure(2);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif

    /** Nodes **/
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(2);
    int node = 0;
    double elementBeginX = 0.;
    double elementBeginY = 0.;

    nodeCoordinates(0) = factorX*nodeCoordinatesFirstElement(0);
    nodeCoordinates(1) = factorY*nodeCoordinatesFirstElement(0);
    //myStructure->NodeCreateDOFs(node, "Displacements", nodeCoordinates);
    myStructure->NodeCreate(node, nodeCoordinates);
    if(PRINTRESULT) std::cout << "create node: " << node << " coordinates: " << nodeCoordinates.at(0,0) <<", "<< nodeCoordinates.at(1,0) << std::endl;
    node++;

    for(int y = 0; y < NumElementsY; y++)
    {
        for (int i = ((y==0)?0:1); i < nodeCoordinatesFirstElement.size(); i++)
        {
            if (node != 1)
            {
                nodeCoordinates(0) = factorX*nodeCoordinatesFirstElement(0);
                nodeCoordinates(1) = factorY*(nodeCoordinatesFirstElement(i) + elementBeginY);
                //myStructure->NodeCreateDOFs(node, "Displacements", nodeCoordinates);
                myStructure->NodeCreate(node, nodeCoordinates);
                if(PRINTRESULT) std::cout << "create node: " << node << " coordinates: " << nodeCoordinates.at(0,0) <<", "<< nodeCoordinates.at(1,0)  << std::endl;
                node++;
            }

            elementBeginX = 0.;
            for(int x = 0; x < NumElementsX; x++)
            {
                for (int j = 1; j < nodeCoordinatesFirstElement.size(); j++)
                {
                    nodeCoordinates(0) = factorX*(nodeCoordinatesFirstElement(j) + elementBeginX);
                    nodeCoordinates(1) = factorY*(nodeCoordinatesFirstElement(i) + elementBeginY);
                    //myStructure->NodeCreateDOFs(node, "Displacements", nodeCoordinates);
                    myStructure->NodeCreate(node, nodeCoordinates);
                    if(PRINTRESULT) std::cout << "create node: " << node << " coordinates: " << nodeCoordinates.at(0,0) <<", "<< nodeCoordinates.at(1,0)  << std::endl;
                    node++;
                }
                elementBeginX += elementSize;
            }
        }
        elementBeginY += elementSize;
    }

    int interpolationType = myStructure->InterpolationTypeCreate("QUAD2D");
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, rElementTypeIdent);
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, rElementTypeIdent);

    /** Elements **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(rNumNodesPerElementInOneDir*rNumNodesPerElementInOneDir);
    int numNodesInRow = NumElementsX*(rNumNodesPerElementInOneDir-1) + 1;
    for(int j = 0; j < NumElementsY; j++)
    {
        for (int i = 0; i < NumElementsX; i++)
        {
            // one element
            for (int k = 0; k < rNumNodesPerElementInOneDir; k++)
                for(int l = 0; l < rNumNodesPerElementInOneDir; l++)
                    elementIncidence(k + l*rNumNodesPerElementInOneDir) = i*(rNumNodesPerElementInOneDir - 1) + k + l*numNodesInRow + j*(rNumNodesPerElementInOneDir-1)*numNodesInRow;

            int myElement = myStructure->ElementCreate(interpolationType, elementIncidence);
            if (PRINTRESULT) std::cout <<  "create element: " << myElement << " nodes: \n" << elementIncidence << std::endl;
        }
    }

    myStructure->Info();
    /** create constitutive law **/
    int myMatLin = myStructure->ConstitutiveLawCreate("LinearElasticEngineeringStress");
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,YoungsModulus);
    myStructure->ConstitutiveLawSetParameterDouble(myMatLin,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,PoissonRatio);

    /** create section **/
    myStructure->ElementTotalConvertToInterpolationType(1.e-6,10);
    //myStructure->InterpolationTypeSetIntegrationType(myInterpolationType, rIntegrationTypeIdent, NuTo::IpData::STATICDATA);
    int mySection = myStructure->SectionCreate("Plane_Stress");
    myStructure->SectionSetThickness(mySection,Thickness);

    /** assign constitutive law **/
    myStructure->ElementTotalSetConstitutiveLaw(myMatLin);
    myStructure->ElementTotalSetSection(mySection);

    // Dirichlet
    NuTo::FullVector<double,Eigen::Dynamic> direction(2);

    if (PRINTRESULT) std::cout <<  "Nodes fixed: "<< std::endl;

    direction(0) = 1; direction(1) = 0;
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    if (PRINTRESULT) std::cout <<  0 << std::endl;
    int nodeId = 0;
    direction(1) = 0;
    for(int j = 0; j < NumElementsY; j++)
    {
        for(int l = 1; l < rNumNodesPerElementInOneDir; l++)
        {
            nodeId = l*numNodesInRow + j*(rNumNodesPerElementInOneDir-1)*numNodesInRow;
            myStructure->ConstraintLinearSetDisplacementNode(nodeId, direction, 0.0);
            if (PRINTRESULT) std::cout <<  nodeId << std::endl;
        }
    }

    direction(1) = 1;
    myStructure->ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    if (PRINTRESULT) std::cout <<  0 << ", "<< direction << std::endl;

    // right boundary
    int groupNumberNodesRight = myStructure->GroupCreate("NODES");
    if (PRINTRESULT) std::cout <<  "Pressure on nodes: "<< std::endl;
    nodeId = (NumElementsX-1)*(rNumNodesPerElementInOneDir-1) + (rNumNodesPerElementInOneDir-1);
    if (PRINTRESULT) std::cout <<  nodeId << std::endl;
    myStructure->GroupAddNode(groupNumberNodesRight, nodeId);

    for(int j = 0; j < NumElementsY; j++)
    {
        for(int l = 1; l < rNumNodesPerElementInOneDir; l++)
        {
            nodeId = (NumElementsX-1)*(rNumNodesPerElementInOneDir-1) +
                    (rNumNodesPerElementInOneDir-1) + l*numNodesInRow + j*(rNumNodesPerElementInOneDir-1)*numNodesInRow;

            if(PRINTRESULT) std::cout << nodeId << std::endl;
            myStructure->GroupAddNode(groupNumberNodesRight, nodeId);
        }
    }

    if (PRINTRESULT) std::cout <<  " and elements: \n";
    int groupNumberElementsRight = myStructure->GroupCreate("ELEMENTS");
    for (int i = 0; i < NumElementsY; i++)
    {
        int elementId = NumElementsX*(i+1)-1;
        myStructure->GroupAddElement(groupNumberElementsRight, elementId);
        if(PRINTRESULT) std::cout << elementId  << std::endl;
    }

    //NuTo::FullVector<double, 2> ForceVecRight({Force,0.});
    myStructure->SetNumLoadCases(1);
    //myStructure->LoadSurfaceConstDirectionCreate2D(0, groupNumberElementsRight, groupNumberNodesRight, ForceVecRight);
    myStructure->LoadSurfacePressureCreate2D(0, groupNumberElementsRight, groupNumberNodesRight, -Stress);

    myStructure->Info();


    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    return myStructure;
}

void solve(NuTo::Structure *myStructure, double solution, double tol = 1.e-6)
{
    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;
    myStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

    stiffnessMatrixCSRVector2.IsSymmetric();

    // build global external load vector
    NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
    myStructure->BuildGlobalExternalLoadVector(0,extForceVector);

    // calculate right hand side
    NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

    #if defined HAVE_PARDISO
    NuTo::SparseDirectSolverPardiso mySolverPardiso(myStructure->GetNumProcessors());
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorPardiso;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverPardiso.Solve(stiffnessMatrix, rhsVector, displacementVectorPardiso);
    // write displacements to node
    myStructure->NodeMergeActiveDofValues(displacementVectorPardiso);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorPardiso;
    myStructure->BuildGlobalGradientInternalPotentialVector(intForceVectorPardiso);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorPardiso = extForceVector - intForceVectorPardiso;
    std::cout << "residual: " << residualVectorPardiso.Norm() << std::endl;
    #elif defined HAVE_MUMPS
    NuTo::SparseDirectSolverMUMPS mySolverMUMPS;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorMUMPS;
    stiffnessMatrix.SetOneBasedIndexing();

    mySolverMUMPS.Solve(stiffnessMatrix, rhsVector, displacementVectorMUMPS);
    // write displacements to node
    myStructure->NodeMergeActiveDofValues(displacementVectorMUMPS);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorMUMPS;
    myStructure->BuildGlobalGradientInternalPotentialVector(intForceVectorMUMPS);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorMUMPS = extForceVector - intForceVectorMUMPS;
    std::cout << "residual: " << residualVectorMUMPS.Norm() << std::endl;
    #else
    std::cout << "Solver not available - can't solve system of equations " << std::endl;
    #endif

    #ifdef ENABLE_VISUALIZE
    // visualize results
    int visualizationGroup = myStructure->GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure->GroupAddElementsTotal(visualizationGroup);

    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    myStructure->ExportVtkDataFileElements("LobattoTruss1D2N.vtk");
    #endif // ENABLE_VISUALIZE

//    double DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    double nodeDisp = myStructure->NodeGetNodePtr(myStructure->GetNumNodes()-1)->GetDisplacement(0);

    std::cout << "Displacement last node FEM:\n" << nodeDisp << std::endl;
    std::cout << "Displacement last node analytical:\n" << solution << std::endl;

    if (fabs(nodeDisp - solution)/fabs(solution) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }
}


int main()
{
    double DisplacementCorrectSerialization1D;
    double DisplacementCorrectSerialization2D;

    // 3 IP
    NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
    {
        NuTo::Structure* myStructure;

        NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;

        for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
        nodeCoordinates += ones;

        // 3Nodes 1D
        myStructure = buildStructure1D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, 10, 2, DisplacementCorrectSerialization1D);
        solve(myStructure, DisplacementCorrectSerialization1D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut1D3NLobatto)**************************\n";
        myStructure->Save("StructureOut1D3NLobatto", "BINARY");
#endif
        // 3Nodes 2D
        myStructure = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, 4, 2, DisplacementCorrectSerialization2D);
        solve(myStructure, DisplacementCorrectSerialization2D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut2D3NLobatto)**************************\n";
        myStructure->Save("StructureOut2D3NLobatto", "XML");
#endif
    }

#ifdef ENABLE_SERIALIZATION
    {
        NuTo::Structure *myStructureImported1D = new NuTo::Structure(1);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut1D3NLobatto)**************************\n";
        myStructureImported1D->Restore("StructureOut1D3NLobatto", "BINARY");
        std::cout << "\n************************** Solve extracted structure (StructureOut1D3NLobatto) **************************\n";
        solve(myStructureImported1D, DisplacementCorrectSerialization1D);

        NuTo::Structure *myStructureImported2D = new NuTo::Structure(2);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut2D3NLobatto)**************************\n";
        myStructureImported2D->Restore("StructureOut2D3NLobatto", "XML");
        std::cout << "\n************************** Solve extracted structure (StructureOut2D3NLobatto) **************************\n";
        solve(myStructureImported2D, DisplacementCorrectSerialization2D);
    }
#endif


    // 4 IP
    ones.resize(4); ones.fill(1);
    nodeCoordinates.resize(4);
    {
        NuTo::Structure* myStructure;

        NuTo::IntegrationType1D2NLobatto4Ip Lobatto1D2N4Ip;
        for (int i = 0; i < 4; i++) Lobatto1D2N4Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
        nodeCoordinates+=ones;


        // 4Nodes 1D
        myStructure = buildStructure1D(NuTo::Interpolation::eTypeOrder::LOBATTO3, 4, nodeCoordinates, 10, 2, DisplacementCorrectSerialization1D);
        solve(myStructure, DisplacementCorrectSerialization1D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut1D4NLobatto)**************************\n";
        myStructure->Save("StructureOut1D4NLobatto", "BINARY");
#endif
        // 3Nodes 2D
        myStructure = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO3, 4, nodeCoordinates, 4, 2, DisplacementCorrectSerialization2D);
        solve(myStructure, DisplacementCorrectSerialization2D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut2D4NLobatto)**************************\n";
        myStructure->Save("StructureOut2D4NLobatto", "XML");
#endif
    }

#ifdef ENABLE_SERIALIZATION
    {
        NuTo::Structure *myStructureImported1D = new NuTo::Structure(1);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut1D4NLobatto)**************************\n";
        myStructureImported1D->Restore("StructureOut1D4NLobatto", "BINARY");
        std::cout << "\n************************** Solve extracted structure (StructureOut1D4NLobatto) **************************\n";
        solve(myStructureImported1D, DisplacementCorrectSerialization1D);

        NuTo::Structure *myStructureImported2D = new NuTo::Structure(2);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut2D4NLobatto)**************************\n";
        myStructureImported2D->Restore("StructureOut2D4NLobatto", "XML");
        std::cout << "\n************************** Solve extracted structure (StructureOut2D4NLobatto) **************************\n";
        solve(myStructureImported2D, DisplacementCorrectSerialization2D);
    }
#endif


    // 5 IP
    ones.resize(5); ones.fill(1);
    nodeCoordinates.resize(5);
    {
        NuTo::Structure* myStructure;

        NuTo::IntegrationType1D2NLobatto5Ip Lobatto1D2N5Ip;
        for (int i = 0; i < 5; i++) Lobatto1D2N5Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
        nodeCoordinates+=ones;


        // 5Nodes 1D
        myStructure = buildStructure1D(NuTo::Interpolation::eTypeOrder::LOBATTO4, 5, nodeCoordinates, 10, 2, DisplacementCorrectSerialization1D);
        solve(myStructure, DisplacementCorrectSerialization1D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut1D5NLobatto)**************************\n";
        myStructure->Save("StructureOut1D5NLobatto", "XML");
#endif
        // 5Nodes 2D
        myStructure = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO4, 5, nodeCoordinates, 4, 2, DisplacementCorrectSerialization2D);
        solve(myStructure, DisplacementCorrectSerialization2D);

#ifdef ENABLE_SERIALIZATION
        std::cout << "\n************************** Saving a NuTo-Structure (StructureOut2D5NLobatto)**************************\n";
        myStructure->Save("StructureOut2D5NLobatto", "BINARY");
#endif
    }

#ifdef ENABLE_SERIALIZATION
    {
        NuTo::Structure *myStructureImported1D = new NuTo::Structure(1);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut1D5NLobatto)**************************\n";
        myStructureImported1D->Restore("StructureOut1D5NLobatto", "XML");
        std::cout << "\n************************** Solve extracted structure (StructureOut1D5NLobatto) **************************\n";
        solve(myStructureImported1D, DisplacementCorrectSerialization1D);

        NuTo::Structure *myStructureImported2D = new NuTo::Structure(2);
        std::cout << "\n************************** Extracting a NuTo-Structure (StructureOut2D5NLobatto)**************************\n";
        myStructureImported2D->Restore("StructureOut2D5NLobatto", "BINARY");
        std::cout << "\n************************** Solve extracted structure (StructureOut2D5NLobatto) **************************\n";
        solve(myStructureImported2D, DisplacementCorrectSerialization2D);
    }
#endif

    return 0;
}
