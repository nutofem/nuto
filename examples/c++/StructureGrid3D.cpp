#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    //   int readFlag = false;

    double PoissonsRatio = 0.2;
    //for local base stiffness matrix
    double YoungsModulus = 1.;

    unsigned int numVoxel;
    const double* voxelSpacing;
    const double* gridOrigin;
    const unsigned int* gridDimension;

    // create structure
    try
    {
        NuTo::StructureGrid myGrid(3);

        // read entries
        myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("/home/fuhlrott/develop/nuto/examples/c++/InputStructureGrid3D");
        numVoxel=myGrid.GetNumVoxels();
        voxelSpacing=myGrid.GetVoxelSpacing();
        gridOrigin=myGrid.GetGridOrigin();
        gridDimension=myGrid.GetGridDimension();
        std::cout<<__FILE__<<"  variab: spac: "<<voxelSpacing[0]<< " gridDim: "<<gridDimension[0]<<std::endl;
        NuTo::FullMatrix<int> imageValues (numVoxel,1);
        imageValues.NuTo::FullMatrix<int>::ImportFromVtkASCIIFile( "/home/fuhlrott/develop/nuto/examples/c++/InputStructureGrid3D");

        std::cout<<"first value "<< imageValues(0,0) << std::endl;
        std::cout<<"numVoxel"<< numVoxel << std::endl;

        //RB
        //double Force = 1.;
        bool EnableDisplacementControl = true;
        double BoundaryDisplacement = 0.1;

        //calculate one element stiffness matrix with E=1
        NuTo::Structure myHelpStruc(3);
        // create material law
        int myMat=myHelpStruc.ConstitutiveLawCreate("LinearElastic");
        myHelpStruc.ConstitutiveLawSetPoissonsRatio(myMat, PoissonsRatio);
        myHelpStruc.ConstitutiveLawSetYoungsModulus(myMat, YoungsModulus);

        // create nodes
        NuTo::FullMatrix<double> nodeCoordinates(3, 1);
        NuTo::FullMatrix<int> elementIncidence(8,1);
        int count = 0;

         // for first voxel

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(0,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(1,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(2,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -voxelSpacing[2] / 2;
        elementIncidence(3,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(4,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(5,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(6,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        nodeCoordinates(0, 0) = -voxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = voxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = voxelSpacing[2] / 2;
        elementIncidence(7,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

        // first element create

        elementIncidence.Info();
        int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
        //myHelpStruc.ElementSetConstitutiveLaw(element,"Material1");
        myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);

        myHelpStruc.NodeInfo(0);

        myHelpStruc.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;
        myHelpStruc.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);

        stiffnessMatrix.RemoveZeroEntries(0,1e-14);
        stiffnessMatrix.Info();
        //NuTo::FullMatrix A(stiffnessMatrix);
        //A.WriteToFile("$HOME/develop/nuto/stiffnessMatrix.txt"," ");

        //grid structure create
        myGrid.CreateNodeGrid("DISPLACEMENTS");

        //myGrid.CreateNodeGrid("COORDINATES");
        std::cout<<"NodeGrid created"<<std::endl;
        int numNodes=myGrid.GetNumNodes();
        std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;
        std::cout<<"knoten "<<myGrid.NodeGetId(numNodes-1)<<std::endl; //Gitterknoten


        //set Modul for each color
        NuTo::FullMatrix<double> myMapColorModul(255,1);
        count=0;
        for(count=0;count<131;count++)
            myMapColorModul(count,0)=0;
        for(count=131;count<161;count++)
            myMapColorModul(count,0)=8300.;
        for (count=161;count<255;count++)
            myMapColorModul(count,0)=11500.;
        //myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");
         myGrid.CreateElementGrid(stiffnessMatrix,myMapColorModul,"VOXEL8N");
        std::cout<<"ElementGrid created"<<std::endl;

        // boundary conditions
        int NumElementsX = (int) gridDimension[0];
        int NumElementsY = (int) gridDimension[1];
        int NumElementsZ = (int) gridDimension[2];

        NuTo::FullMatrix<double> direction(3,1);
        direction(0,0)= 1;
        direction(1,0)= 0;
        direction(2,0)= 0;
        std::cout<<"num nodes "<<myGrid.GetNumNodes()<<std::endl;

      for(int zCount = 0; zCount < NumElementsZ + 1; zCount++)
        {
            for(int yCount = 0; yCount < NumElementsY + 1; yCount++)
            {
                int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1);
                int myNodeNumber;
                int flag=0;
                try
                {
                    myNodeNumber=myGrid.NodeGetNodeNumberFromId(node); //node from type NodeGridCoordinates
                }
                catch(NuTo::MechanicsException& e)
                {
                    flag=1;
                }
                if(flag==0)
                     myGrid.ConstraintSetDisplacementNode(myNodeNumber, direction, 0.0);
            }
        }
        direction(0,0)= 0;
        direction(1,0)= 0;
        direction(2,0)= 1;
        myGrid.ConstraintSetDisplacementNode(0, direction, 0.0);
        myGrid.ConstraintSetDisplacementNode(NumElementsY * (NumElementsX + 1), direction, 0.0);
        direction(0,0)= 0;
        direction(1,0)= 1;
        direction(2,0)= 0;
        myGrid.ConstraintSetDisplacementNode(0, direction, 0.0);

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
                    int myNodeNumber;
                    int flag=0;
                    try
                    {
                        myNodeNumber=myGrid.NodeGetNodeNumberFromId(node); //node from type NodeGridCoordinates
                    }
                    catch(NuTo::MechanicsException& e)
                    {
                        flag=1;
                    }
                   if(flag==0)
                        myGrid.ConstraintSetDisplacementNode(myNodeNumber, direction, BoundaryDisplacement);
                }
            }
         }
        else
        {
            std::cout << "Load control" << "not implemented"<<std::endl;
         /*   //! @TODO: Add special configurations if edge node not exist
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
                int myNodeNumber;
                int flag=0;
                try
                {
                  myNodeNumber=myGrid.NodeGetNodeNumberFromId(node); //node from type NodeGridCoordinates
                }
                catch(NuTo::MechanicsException& e)
                {
                    std::cout<<__FILE__<<__LINE__<<"node existiert nicht"<<node<<std::endl;
                    flag=1;
                }
                if(flag==0)
                {
                    std::cout<<__FILE__<<__LINE__<<"myNodeNumber"<<myNodeNumber<<"  von  "<<myGrid.GetNumNodes()<<std::endl;
                    myGrid.LoadCreateNodeForce(myNodeNumber, direction, nodeForce);
                }
                for(int yCount = 1; yCount < NumElementsY; yCount++)
                {
                    node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
                    std::cout << "apply force to node: " << node << " force: " << 2 * nodeForce << std::endl;

                    myGrid.LoadCreateNodeForce(node, direction, 2 * nodeForce);
                }
                node = (zCount + 1) * (NumElementsX + 1) * (NumElementsY + 1) - 1;
                //std::cout << "apply force to node: " << node << " force: " << nodeForce << std::endl;
                myGrid.LoadCreateNodeForce(node, direction, nodeForce);
            }
            */

        }
        // start analysis
        std::cout<<__FILE__<<" "<<__LINE__<<"  start analysis"<<std::endl;
        // build global dof numbering
        myGrid.NodeBuildGlobalDofs();
        std::cout<<__FILE__<<" "<<__LINE__<<"  glob dofs"<<myGrid.NodeGetNumberGlobalDofs()<<std::endl;
        std::cout<<__FILE__<<" "<<__LINE__<<" active dofs"<<myGrid.NodeGetNumberActiveDofs()<<std::endl;
/*
        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> globStiffnessMatrix;
        myGrid.BuildGlobalCoefficientMatrix0(globStiffnessMatrix, dispForceVector);
        globStiffnessMatrix.RemoveZeroEntries(0,1e-14);
        //NuTo::FullMatrix<double> A(stiffnessMatrix);
        //A.WriteToFile("stiffnessMatrix.txt"," ");
        //stiffnessMatrix.Info();
        dispForceVector.Info();

        // build global external load vector
        NuTo::FullMatrix<double> extForceVector;
        myGrid.BuildGlobalExternalLoadVector(extForceVector);
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
        myGrid.NodeMergeActiveDofValues(displacementVector);

        // calculate residual
        NuTo::FullMatrix<double> intForceVector;
        myGrid.BuildGlobalGradientInternalPotentialVector(intForceVector);
        NuTo::FullMatrix<double> residualVector = extForceVector - intForceVector;
        std::cout << "residual: " << residualVector.Norm() << std::endl;
*/
        // visualize results
//        myGrid.ExportVtkDataFile("StrukturedGrid3D.vtk","DISPLACEMENTS");


    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureGrid3D] error.");
        std::cerr<<e.ErrorMessage()<<std::endl;
        return -1;
    }
    catch(...)
    {
        throw NuTo::MechanicsException("[NuTo::StructureGrid3D] Unknown error catched .");
    }
    return 0;

}
