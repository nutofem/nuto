#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

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

#define PRINTRESULT false

#include <string>

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
                     int rNumNodesPerElementInOneDir,
                     NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                     int NumElementsX, int NumElementsY)
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

    /** Structure 2D **/
    NuTo::Structure *myStructure = new NuTo::Structure(2);

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

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    myStructure->Info();

    return myStructure;
}

void solve2d(NuTo::Structure* myStructure,
             int NumElementsX,
             int NumElementsY,
             double  YoungsModulus = 20000.,
             double  Stress = 10.,
             double  Length = 10.0,
             double  tol = 1.e-6)
{
    /** start analysis **/
    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values

    //NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2(0,0);
    NuTo::SparseMatrixCSRVector2Symmetric<double> stiffnessMatrixCSRVector2Symmetric(0,0);

    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector;

    //myStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    myStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2Symmetric, dispForceVector);

    //NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);
    NuTo::SparseMatrixCSRSymmetric<double> stiffnessMatrixSymmetric (stiffnessMatrixCSRVector2Symmetric);

    // build global external load vector
    NuTo::FullVector<double,Eigen::Dynamic> extForceVector;
    myStructure->BuildGlobalExternalLoadVector(0,extForceVector);

    // calculate right hand side
    NuTo::FullVector<double,Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

    //stiffnessMatrix.SetOneBasedIndexing();
    stiffnessMatrixSymmetric.SetOneBasedIndexing();


#ifdef HAVE_PARDISO
    NuTo::SparseDirectSolverPardiso mySolverPardiso(myStructure->GetNumProcessors());
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorPardiso;
    //mySolverPardiso.Solve(stiffnessMatrix, rhsVector, displacementVectorPardiso);
    mySolverPardiso.Solve(stiffnessMatrixSymmetric, rhsVector, displacementVectorPardiso);
    // write displacements to node
    myStructure->NodeMergeActiveDofValues(displacementVectorPardiso);
    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorPardiso;
    myStructure->BuildGlobalGradientInternalPotentialVector(intForceVectorPardiso);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorPardiso = extForceVector - intForceVectorPardiso;
    std::cout << "residual: " << residualVectorPardiso.Norm() << std::endl;
    // visualize results
//    myStructure->AddVisualizationComponentDisplacements();
//    myStructure->AddVisualizationComponentEngineeringStrain();
//    myStructure->AddVisualizationComponentEngineeringStress();
    myStructure->ExportVtkDataFileElements("Lobatto2D.vtk");
#elif defined HAVE_MUMPS
    NuTo::SparseDirectSolverMUMPS mySolverMUMPS;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVectorMUMPS;

    mySolverMUMPS.Solve(stiffnessMatrixSymmetric, rhsVector, displacementVectorMUMPS);
    // write displacements to node
    myStructure->NodeMergeActiveDofValues(displacementVectorMUMPS);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVectorMUMPS;
    myStructure->BuildGlobalGradientInternalPotentialVector(intForceVectorMUMPS);
    NuTo::FullVector<double,Eigen::Dynamic> residualVectorMUMPS = extForceVector - intForceVectorMUMPS;
    std::cout << "residual: " << residualVectorMUMPS.Norm() << std::endl;
#else
    std::cout << "PARDISO ans MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_PARDISO

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStress(6,1);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStressCorrect(6,1);
    EngineeringStressCorrect(0,0) = Stress;

    int numNodes = myStructure->GetNumNodes();

    double DisplacementCorrect = (Stress*Length)/YoungsModulus;
    double nodeDisp = myStructure->NodeGetNodePtr(numNodes-1)->GetDisplacement(0);

    std::cout << "Displacement analytical: " << DisplacementCorrect  << std::endl;
    std::cout << "Displacement FEM: " <<  nodeDisp << std::endl;

    for(int element = 0; element < NumElementsX*NumElementsY; element++)
    {
        if (PRINTRESULT)
        {
            myStructure->ElementGetEngineeringStress(element, EngineeringStress);
            std::cout << "EngineeringStressCorrect" << std::endl;
            EngineeringStressCorrect.Info();
            std::cout << "EngineeringStress" << std::endl;
            EngineeringStress.Info();
        }
    }

    double resultTol = fabs(nodeDisp - DisplacementCorrect)/fabs(DisplacementCorrect);

    if (fabs(resultTol) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration2D] : displacement is not correct");
    }
}


/*
 |>*----*----*----*----*  -->F
*/
NuTo::Structure* buildStructure1D(NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                                  int rNumNodesPerElement,
                                  NuTo::FullVector<double, Eigen::Dynamic>& nodeCoordinatesFirstElement,
                                  int     NumElements = 10,
                                  double  YoungsModulus = 20000.,
                                  double  Area = 100.0*0.1,
                                  double  Length = 1000.0,
                                  double  Force = 1.,
                                  double  tol = 1.e-6)
{
    /** Structure 1D **/
    NuTo::Structure *myStructure = new NuTo::Structure(1);

#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif // OPENMP

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

    int groupNumberNodes = myStructure->GroupCreate("NODES");
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
            myStructure->GroupAddNode(groupNumberNodes, node);
            node++;
        }
        elementBegin+=elementLength;
    }

    int interpolationType = myStructure->InterpolationTypeCreate("TRUSS1D");
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::COORDINATES, rElementTypeIdent);
    myStructure->InterpolationTypeAdd(interpolationType, NuTo::Node::DISPLACEMENTS, rElementTypeIdent);

    int numNodes = node;

    /** create elements **/
    int groupNumberElements = myStructure->GroupCreate("ELEMENTS");
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

        int myElement =  myStructure->ElementCreate(interpolationType, elementIncidence);
        myStructure->GroupAddElement(groupNumberElements, myElement);
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

    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

    return myStructure;
}

void solve1d(NuTo::Structure* myStructure,
           double  YoungsModulus = 20000.,
           double  Area = 100.0*0.1,
           double  Length = 1000.0,
           double  Force = 1.,
           double  tol = 1.e-6)
{
    /** start analysis **/
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

#ifdef HAVE_PARDISO
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

    double DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    int numNodes = myStructure->GetNumNodes();
    double nodeDisp = myStructure->NodeGetNodePtr(numNodes-1)->GetDisplacement(0);

    std::cout << "Displacement last node FEM:\n" << nodeDisp << std::endl;
    std::cout << "Displacement last node analytical:\n" << DisplacementCorrect << std::endl;

    if (fabs(nodeDisp - DisplacementCorrect)/fabs(DisplacementCorrect) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }

}

#ifdef ENABLE_SERIALIZATION
void serialize1d()
{
    {
        std::cout << "*** Serializing a NuTo-Structure to ***\n";
        NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
        NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
        NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;
        for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
        nodeCoordinates += ones;
        NuTo::Structure *myStructure = buildStructure1D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, 2);

//        std::ofstream outFileStream("StructureOut");
//        boost::archive::xml_oarchive outArchivexml(outFileStream);
//        myStructure->save(outArchivexml, 1);
//        outFileStream.close();
        myStructure->Save("StructureOut", "XML");
    }



    {
        std::cout << "\n************************** Extracting a NuTo-Structure from StructureOut **************************\n";
        NuTo::Structure *myStructureImported = new NuTo::Structure(1);
    //    std::ifstream iFileStreamxml("StructureOut");
    //    boost::archive::xml_iarchive inArchivexml(iFileStreamxml);
    //    myStructureImported->load(inArchivexml, 1);
        myStructureImported->Restore("StructureOut", "XML");
        myStructureImported->GroupInfo(3);
        std::cout << "\n\n\n\n\n*** Solve extracted structure ***\n\n\n\n\n";
        solve1d(myStructureImported);
    }
}

void serialize2d()
{
    int numElementsX = 10;
    int numElementsY = 10;
    {
        //  3Nodes 2D
        NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
        NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
        NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;
        ones.resize(3); ones.fill(1);
        nodeCoordinates.resize(3);
        for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
        nodeCoordinates += ones;
        NuTo::Structure* myStructure2D = buildStructure2D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, numElementsX, numElementsY);

//        std::ofstream outFileStream("StructureOut2D");
//        boost::archive::xml_oarchive outArchivexml(outFileStream);
//        myStructure2D->save(outArchivexml, 1);
//        outFileStream.close();
        myStructure2D->Save("StructureOut2D", "XML");
    }

    {
        std::cout << "\n************************** Extracting a NuTo-Structure from StructureOut2D **************************\n";
        NuTo::Structure *myStructureImported = new NuTo::Structure(2);
//        std::ifstream iFileStreamxml("StructureOut2D");
//        boost::archive::xml_iarchive inArchivexml(iFileStreamxml);
//        myStructureImported->load(inArchivexml, 1);
        myStructureImported->Restore("StructureOut2D", "XML");

        std::cout << "\n************************** Solve extracted structure 2D **************************\n";
        solve2d(myStructureImported,numElementsX,numElementsY);
    }
}
#endif

int main(int argc, char* argv[])
{

#ifdef ENABLE_SERIALIZATION
    serialize1d();
    serialize2d();
#endif

    return EXIT_SUCCESS;
}


