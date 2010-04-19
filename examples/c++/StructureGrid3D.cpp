#include <iostream>
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

int main()
{
    //   int readFlag = false;

    double PoissonsRatio = 0.2;
    //for local base stiffness matrix
    double YoungsModulus = 1.;

    unsigned int numVoxel;
    const double* mVoxelSpacing;
    const double* gridOrigin;
    const unsigned int* gridDimension;

    // create structure
    try
    {
        NuTo::StructureGrid myGrid(3);
        std::cout<<"test2 \n";


        // read entries
        myGrid.NuTo::StructureGrid::ImportFromVtkASCIIFileHeader("/home/fuhlrott/develop/nuto/intNutoInputsmall");
        numVoxel=myGrid.GetNumVoxels();
        mVoxelSpacing=myGrid.GetVoxelSpacing();
        gridOrigin=myGrid.GetGridOrigin();
        gridDimension=myGrid.GetGridDimension();
        NuTo::FullMatrix<int> imageValues (numVoxel,1);
        imageValues.NuTo::FullMatrix<int>::ImportFromVtkASCIIFile( "/home/fuhlrott/develop/nuto/intNutoInputsmall");

        std::cout<<"first value "<< imageValues(0,0) << std::endl;
        std::cout<<"numVoxel"<< numVoxel << std::endl;


        //calculate one element stiffness matrix with E=1
        NuTo::Structure myHelpStruc(3);
        // create material law
        myHelpStruc.ConstitutiveLawCreate("Material1", "LinearElastic");
        myHelpStruc.ConstitutiveLawSetPoissonsRatio("Material1", PoissonsRatio);
        myHelpStruc.ConstitutiveLawSetYoungsModulus("Material1", YoungsModulus);

        // create nodes
        NuTo::FullMatrix<double> nodeCoordinates(3, 1);
        int node = 0;
        int count = 0;

        std::cout<<"test4 \n";

         // for first voxel

        nodeCoordinates(0, 0) = -mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = -mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = -mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = -mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = -mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        node++;
        nodeCoordinates(0, 0) = -mVoxelSpacing[0] / 2;
        nodeCoordinates(1, 0) = mVoxelSpacing[1] / 2;
        nodeCoordinates(2, 0) = mVoxelSpacing[2] / 2;
        myHelpStruc.NodeCreate(node, "displacements", nodeCoordinates);
        myHelpStruc.NodeInfo(0);

        // first element create
        // create elements
        // create elements
        NuTo::FullMatrix<int> elementIncidence(8,1);
        int element = 0;
        elementIncidence(0,0) = 0;
        elementIncidence(1,0) = 1;
        elementIncidence(2,0) = 2;
        elementIncidence(3,0) = 3;
        elementIncidence(4,0) = 4;
        elementIncidence(5,0) = 5;
        elementIncidence(6,0) = 6;
        elementIncidence(7,0) = 7;
        std::cout << "element: " << element << " incidence: " << std::endl;
        //elementIncidence.Info();
        myHelpStruc.ElementCreate(element, "Brick8N", elementIncidence);
        myHelpStruc.ElementSetConstitutiveLaw(element,"Material1");

        myHelpStruc.NodeInfo(0);

        // calculate stiffness matrix
        // build global dof numbering
        //myHelpStruc.NodeBuildGlobalDofs();

        // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix;
        NuTo::FullMatrix<double> dispForceVector;

        myHelpStruc.BuildGlobalCoefficientMatrix0(stiffnessMatrix, dispForceVector);
        stiffnessMatrix.RemoveZeroEntries(0,1e-14);
        NuTo::FullMatrix<double> A(stiffnessMatrix);
        A.WriteToFile("$HOME/develop/nuto/stiffnessMatrix.txt"," ");
        stiffnessMatrix.Info();
        //myHelpStuc löschen

        //grid structure create
        myGrid.CreateNodeGrid("DISPLACEMENTS");
        std::cout<<"NodeGrid created"<<std::endl;

        //set Modul for each color
        NuTo::FullMatrix<double> myMapColorModul(255,1);
        count=0;
        for(count=0;count<131;count++)
            myMapColorModul(count,0)=0;
        for(count=131;count<161;count++)
            myMapColorModul(count,0)=8300.;
        for (count=161;count<255;count++)
            myMapColorModul(count,0)=11500.;
        myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");
        myGrid.CreateElementGrid(myMapColorModul,"BRICK8N");
        std::cout<<"ElementGrid created"<<std::endl;
        std::cout<<"knoten "<<myGrid.GetNumNodes()<<std::endl;
        std::cout<<"knoten "<<myGrid.NodeGetID(GetNumNodes()-1)<<std::endl;
    }
    catch (NuTo::MechanicsException& e)
    {
        e.AddMessage("[NuTo::StructureGrid3D) error.");
        std::cout<<e.ErrorMessage()<<std::endl;
        exit;
    }
}
