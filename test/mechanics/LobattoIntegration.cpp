#include <iostream>
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"

#define PRINTRESULT true

/*
 |>*----*----*----*----*  -->F
*/

int buildStructure1D(const std::string& rIntegrationTypeIdent)
{
    /** paramaters **/
    double  YoungsModulus = 20000.;
    double  Area = 100. * 0.1;
    double  Length = 1000.;
    int     NumElements = 10;
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
    NuTo::FullVector<double,Eigen::Dynamic> nodeCoordinates(1);
    for(int node = 0; node < NumElements + 1; node++)
    {
        std::cout << "create node: " << node << " coordinates: " << node * Length/NumElements << std::endl;

        nodeCoordinates(0) = node * Length/NumElements;
        myStructure.NodeCreate(node, "displacements", nodeCoordinates);
    }

    /** create elements **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(2);
    for(int element = 0; element < NumElements; element++)
    {
        std::cout <<  "create element: " << element << " nodes: " << element << "," << element+1 << std::endl;

        elementIncidence(0) = element;
        elementIncidence(1) = element + 1;
        int myElement = myStructure.ElementCreate("Truss1D2N", elementIncidence);
        myStructure.ElementSetSection(myElement, Section);
        myStructure.ElementSetConstitutiveLaw(myElement, Material);
        myStructure.ElementSetIntegrationType(myElement, rIntegrationTypeIdent, "NOIPDATA");
    }

    /** set boundary conditions and loads **/
    NuTo::FullVector<double,Eigen::Dynamic> direction(1);
    direction(0) = 1;
    // first node is fixed
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    myStructure.LoadCreateNodeForce(1, NumElements, direction, Force);

    /** start analysis **/
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

    // solve
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
    stiffnessMatrix.SetOneBasedIndexing();

#ifdef HAVE_MUMPS
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVector);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
    std::cout << "residual: " << residualVector.Norm() << std::endl;

    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFile("LobattoTruss1D2N.vtk");

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS

    double DisplacementCorrect;
    DisplacementCorrect = (Force*Length)/(Area*YoungsModulus);

    if(PRINTRESULT)
    {
        std::cout << "Displacement last node FEM:\n" << displacementVector(NumElements-1) << std::endl;
        std::cout << "Displacement last node analytical:\n" << DisplacementCorrect << std::endl;
    }

    if (fabs(DisplacementCorrect - displacementVector(NumElements-1))/fabs(DisplacementCorrect) > tol)
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
int buildStructure2D(const std::string& rIntegrationTypeIdent, int numIP)
{
    /** parameters **/
    double  YoungsModulus = 20000.;
    double  PoissonRatio = 0.3;
    double  Height = 100.;
    // there should be no dependency, because the pressure at the boundary is predefined
    double  Thickness = 2.123548;
    double  Length = 1000.;
    int     NumElements = 1;
    double  Stress = 10.;
    double  tol = 1.e-6;

    /** Structure 2D **/
    NuTo::Structure myStructure(2);

    /** create nodes **/
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> nodeCoordinates(2, 2*(NumElements+1));

    /** nodes **/
    for (int column = 0; column < 2*(NumElements+1); column+=2)
    {
        nodeCoordinates(0,column) = column*Length/NumElements/2.;
        nodeCoordinates(1,column) = Height;

        nodeCoordinates(0,column+1) = column*Length/NumElements/2.;
        nodeCoordinates(1,column+1) = 0.;
    }
    myStructure.NodesCreate("displacements",nodeCoordinates);

    /** nodes numbers **/
    NuTo::FullVector<int,Eigen::Dynamic> elementIncidence(4);
    for (int column = 0; column < NumElements; column++)
    {
        elementIncidence(0) = 2*column;
        elementIncidence(1) = 2*column+1;
        elementIncidence(2) = 2*column+3;
        elementIncidence(3) = 2*column+2;

        int myElement = myStructure.ElementCreate("PLANE2D4N", elementIncidence);
        myStructure.ElementSetIntegrationType(myElement, rIntegrationTypeIdent, "NOIPDATA");
    }

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
    direction(0) = 1;
    direction(1) = 0;
    myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0);
    myStructure.ConstraintLinearSetDisplacementNode(0, direction, 0.0);
    // fix node 1 in y
    direction(0) = 0;
    direction(1) = 1;
    myStructure.ConstraintLinearSetDisplacementNode(1, direction, 0.0);

    // right boundary
    int groupNumberNodesRight = myStructure.GroupCreate("NODES");
    myStructure.GroupAddNode(groupNumberNodesRight, 2*NumElements+1);
    myStructure.GroupAddNode(groupNumberNodesRight, 2*NumElements);

    int groupNumberElementsRight = myStructure.GroupCreate("ELEMENTS");
    myStructure.GroupAddElement(groupNumberElementsRight, NumElements-1);

    //NuTo::FullVector<double, 2> ForceVecRight({Force,0.});
    //myStructure.LoadSurfaceConstDirectionCreate2D(1, groupNumberElementsRight, groupNumberNodesRight, ForceVecRight);
    myStructure.LoadSurfacePressureCreate2D(1, groupNumberElementsRight, groupNumberNodesRight, -Stress);

    myStructure.Info();

    /** start analysis **/
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

    // solve
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullVector<double,Eigen::Dynamic> displacementVector;
    stiffnessMatrix.SetOneBasedIndexing();

#ifdef HAVE_MUMPS
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVector);

    // calculate residual
    NuTo::FullVector<double,Eigen::Dynamic> intForceVector;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    NuTo::FullVector<double,Eigen::Dynamic> residualVector = extForceVector - intForceVector;
    std::cout << "residual: " << residualVector.Norm() << std::endl;

    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.ExportVtkDataFileElements("Lobatto2D4N.vtk");

#else
    std::cout << "MUMPS not available - can't solve system of equations " << std::endl;
#endif // HAVE_MUMPS

    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStress(6,1);
    NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> EngineeringStressCorrect(6,1);
    EngineeringStressCorrect(0,0) = Stress;

    double displacementX = 0.5;

    if (PRINTRESULT)
    {
        std::cout << "Displacement analytical: " << displacementX << std::endl;
        std::cout << "Displacement FEM: " << displacementVector(2*2*(NumElements+1) - 5) << std::endl;
    }

    if ( fabs(displacementX - displacementVector(2*2*(NumElements+1) - 5))/fabs(displacementX) > tol )
    {
        throw NuTo::Exception("[LobattoIntegration(" + rIntegrationTypeIdent + ")] : Displacement is not correct.");
    }

    for(int element = 0; element < NumElements; element++)
    {
        myStructure.ElementGetEngineeringStress(element, EngineeringStress);
        if (PRINTRESULT)
        {
            std::cout << "Displacement" << displacementVector << std::endl;
            std::cout << "EngineeringStressCorrect" << std::endl;
            EngineeringStressCorrect.Info();
            std::cout << "EngineeringStress" << std::endl;
            EngineeringStress.Info();
        }

        for (int countIP=0; countIP<numIP; countIP++)
        {
            if ((EngineeringStress.col(countIP) - EngineeringStressCorrect).cwiseAbs().maxCoeff() > tol)
            {
                throw NuTo::Exception("[LobattoIntegration(" + rIntegrationTypeIdent + ")] : Stress is not correct.");
            }
        }
    }

    return 0;
}

int main()
{
      if (!buildStructure2D("2D4NLobatto9Ip",9) && !buildStructure2D("2D4NLobatto16Ip",16) && !buildStructure2D("2D4NLobatto25Ip",25) &&
          !buildStructure1D("1D2NLobatto3Ip") && !buildStructure1D("1D2NLobatto4Ip") && !buildStructure1D("1D2NLobatto5Ip"))
        return 0;
    else
        return -1;
    return 0;
}
