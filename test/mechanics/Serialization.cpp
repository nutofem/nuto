#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/nodes/NodeDof.h"
//#include "nuto/mechanics/nodes/NodeDof_Def.h"


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
 |>*----*----*----*----*  -->F
*/

NuTo::Structure* buildStructure1D(NuTo::Interpolation::eTypeOrder rElementTypeIdent,
                     //NuTo::IntegrationType::eIntegrationType rIntegrationTypeIdent,
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
    NuTo::Structure *myStructure = new NuTo::Structure(1);

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

            {
                std::string file = "NodeOut.peter";
                std::cout << "\n\n\nWriting a node to " << file << "\n\n\n";
                const NuTo::NodeBase* myNodePtr = myStructure->NodeGetNodePtr(node);
                // serialize it
                std::ofstream outFileStream(file);
                boost::archive::text_oarchive outArchive(outFileStream);

                NuTo::NodeDof<1, 0, 0, 0, 0, 0, 0, 0, 0, 0> myNodeDofDef;
                const NuTo::NodeBase* myNodePtrLoc = &myNodeDofDef;

                outArchive << const_cast<NuTo::NodeBase*&>(myNodePtrLoc);

                outArchive << const_cast<NuTo::NodeBase*&>(myNodePtr);
            }

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

#ifdef _OPENMP
    int numThreads = 4;
    myStructure->SetNumProcessors(numThreads);
#endif // OPENMP

#if defined HAVE_PARDISO

#ifndef _OPENMP
    int numThreads = 1;
#endif

    NuTo::SparseDirectSolverPardiso mySolverPardiso(numThreads);
#endif // PARDISO
    myStructure->CalculateMaximumIndependentSets();
    myStructure->NodeBuildGlobalDofs();

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
    int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure->GroupAddElementsTotal(visualizationGroup);

    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure->AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);

    myStructure->ExportVtkDataFileElements("LobattoTruss1D2N.vtk");
#endif // ENABLE_VISUALIZE

    double DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    double nodeDisp = myStructure->NodeGetNodePtr(numNodes-1)->GetDisplacement(0);

    std::cout << "Displacement last node FEM:\n" << nodeDisp << std::endl;
    std::cout << "Displacement last node analytical:\n" << DisplacementCorrect << std::endl;

    if (fabs(nodeDisp - DisplacementCorrect)/fabs(DisplacementCorrect) > tol)
    {
        throw NuTo::Exception("[LobattoIntegration] : displacement is not correct");
    }

    return myStructure;
}

int main(int argc, char* argv[])
{
    std::cout << "*** Serializing a NuTo-Structure to ***\n";

    NuTo::FullVector<double, Eigen::Dynamic> ones(3); ones.fill(1);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(3);
    NuTo::IntegrationType1D2NLobatto3Ip Lobatto1D2N3Ip;
    for (int i = 0; i < 3; i++) Lobatto1D2N3Ip.GetLocalIntegrationPointCoordinates1D(i, nodeCoordinates(i));
    nodeCoordinates += ones;

    NuTo::Structure *myStructure = buildStructure1D(NuTo::Interpolation::eTypeOrder::LOBATTO2, 3, nodeCoordinates, 10);

    std::string file = "StructureOut";
    std::ofstream outFileStream(file);
    boost::archive::text_oarchive outArchive(outFileStream);
    myStructure->save(outArchive, 1);

    myStructure->Save("StructureOut", "TEXT");

    std::cout << "*** Extracting a NuTo-Structure from StructureOut ***\n";

    NuTo::Structure* myStructureImported;

    myStructureImported->Restore("StructureOut", "TEXT");

    return EXIT_FAILURE;
}


