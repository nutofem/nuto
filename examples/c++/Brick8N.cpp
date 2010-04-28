// $Id$

#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    int readFlag = false;
    if(readFlag)
    {
    NuTo::FullMatrix<double> a;
    a.ReadFromFile("stiffnessMatrix.txt");
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(a);
    NuTo::FullMatrix<double> rhsVector;
    rhsVector.ReadFromFile("rhsVector.txt");
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullMatrix<double> displacementVector;
    //stiffnessMatrix.SetOneBasedIndexing();
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
    displacementVector.WriteToFile("disp.txt"," ");
    a = rhsVector - stiffnessMatrix * displacementVector;
    std::cout << "residual: " << a.Norm() << std::endl;
    }
    else
    {
        // definitions
        double YoungsModulus = 20000.;
        double PoissonsRatio = 0.2;
        double Width = 1000.;
        double Height = 1000.;
        double Length = 1000.;
        int NumElementsX = 6;
        int NumElementsY = 6;
        int NumElementsZ = 6;
        double Force = 1.;
        bool EnableDisplacementControl = true;
        double BoundaryDisplacement = 0.1;

        // create structure
        NuTo::Structure myStructure(3);

        // create material law
        int Material1 = myStructure.ConstitutiveLawCreate("LinearElastic");
        myStructure.ConstitutiveLawSetYoungsModulus(Material1, YoungsModulus);
        myStructure.ConstitutiveLawSetPoissonsRatio(Material1, PoissonsRatio);

        // create nodes
        NuTo::FullMatrix<double> nodeCoordinates(3,1);
        int node = 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            nodeCoordinates(2,0) = (double)zCount * Height/(double)NumElementsZ;
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                nodeCoordinates(1,0) = (double)yCount * Width/(double)NumElementsY;
                for(int xCount = 0; xCount < NumElementsX + 1; xCount++)
                {
                    nodeCoordinates(0,0) = (double)xCount * Length/(double)NumElementsX;
                    //std::cout << "node: " << node << " coordinates: " << nodeCoordinates.GetValue(0,0) << "," << nodeCoordinates.GetValue(1,0) << "," << nodeCoordinates.GetValue(2,0) << std::endl;
                    myStructure.NodeCreate(node, "displacements", nodeCoordinates);
                    node ++;
                }
            }
        }

        // create elements
        NuTo::FullMatrix<int> elementIncidence(8,1);
        int element = 0;
        for(int zCount = 0; zCount < NumElementsZ; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY; yCount++)
            {
                for(int xCount = 0; xCount < NumElementsX; xCount++)
                {
                    int node1 = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + xCount;
                    elementIncidence(0,0) = node1;
                    elementIncidence(1,0) = node1 + 1;
                    elementIncidence(2,0) = node1 + NumElementsX + 2;
                    elementIncidence(3,0) = node1 + NumElementsX + 1;
                    elementIncidence(4,0) = node1 + (NumElementsX + 1) * (NumElementsY + 1);
                    elementIncidence(5,0) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + 1;
                    elementIncidence(6,0) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 2;
                    elementIncidence(7,0) = node1 + (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX + 1;
                    //std::cout << "element: " << element << " incidence: " << std::endl;
                    //elementIncidence.Info();
                    myStructure.ElementCreate(element, "Brick8N", elementIncidence);
                    myStructure.ElementSetConstitutiveLaw(element,Material1);
                    element ++;
                }
            }
        }

        // boundary conditions
        NuTo::FullMatrix<double> direction(3,1);
        direction(0,0)= 1;
        direction(1,0)= 0;
        direction(2,0)= 0;
        for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1);
                //std::cout << "node: " << node << std::endl;
                myStructure.ConstraintSetDisplacementNode(node, direction, 0.0);
            }
        }
        direction(0,0)= 0;
        direction(1,0)= 0;
        direction(2,0)= 1;
        myStructure.ConstraintSetDisplacementNode(0, direction, 0.0);
        myStructure.ConstraintSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0);
        direction(0,0)= 0;
        direction(1,0)= 1;
        direction(2,0)= 0;
        myStructure.ConstraintSetDisplacementNode(0, direction, 0.0);

        // apply nodes
        if(EnableDisplacementControl)
        {
            std::cout << "Displacement control" << std::endl;
            // boundary displacments
            direction(0,0)= 1;
            direction(1,0)= 0;
            direction(2,0)= 0;
            for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
            {
                //std::cout << zCount << std::endl;
                for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
                {
                    int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    //std::cout << "node: " << node << std::endl;
                    myStructure.ConstraintSetDisplacementNode(node, direction, BoundaryDisplacement);
                }
            }
        }
        else
        {
            std::cout << "Load control" << std::endl;
            // apply load to nodes
            direction(0,0)= 1;
            direction(1,0)= 0;
            direction(2,0)= 0;
            for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
            {
                double nodeForce;
                if(zCount == 0 || zCount == NumElementsZ)
                {
                    nodeForce = Force / (4 *NumElementsY * NumElementsZ);
                }
                else
                {
                    nodeForce = Force / (2 *NumElementsY * NumElementsZ);
                }
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + NumElementsX;
                //std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
                myStructure.LoadCreateNodeForce(node, direction, nodeForce);
                for(int yCount = 1; yCount < NumElementsY; yCount++)
                {
                    node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    //std::cout << "apply force to node: " << node << " force: " << 2 * nodeForce << std::endl;
                    myStructure.LoadCreateNodeForce(node, direction, 2 * nodeForce);
                }
                node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1;
                //std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
                myStructure.LoadCreateNodeForce(node, direction, nodeForce);
            }
        }

        // start analysis
        // build global dof numbering
        myStructure.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;
        myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
        stiffnessMatrix.RemoveZeroEntries(0,1e-14);
        //NuTo::FullMatrix<double> A(stiffnessMatrix);
        //A.WriteToFile("stiffnessMatrix.txt"," ");
        //stiffnessMatrix.Info();
        //dispForceVector.Info();

        // build global external load vector
        NuTo::FullMatrix<double> extForceVector;
        myStructure.BuildGlobalExternalLoadVector(extForceVector);
        //extForceVector.Info();

        // calculate right hand side
        NuTo::FullMatrix<double> rhsVector = dispForceVector + extForceVector;
        rhsVector.WriteToFile("rhsVector.txt"," ");

        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        NuTo::FullMatrix<double> displacementVector;
        stiffnessMatrix.SetOneBasedIndexing();
        mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);
        displacementVector.WriteToFile("displacementVector.txt"," ");

        // write displacements to node
        myStructure.NodeMergeActiveDofValues(displacementVector);

        // calculate residual
        NuTo::FullMatrix<double> intForceVector;
        myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
        NuTo::FullMatrix<double> residualVector = extForceVector - intForceVector;
        std::cout << "residual: " << residualVector.Norm() << std::endl;

        // visualize results
        myStructure.ExportVtkDataFile("Brick8N.vtk","displacements engineering_strain engineering_stress");
    }
    return 0;
}

