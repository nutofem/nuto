// $Id$

//! @author Andrea Keßler, ISM
//! @date December 2011
//! @brief ... example for conjugate gradient grid routine with reduced vector size to existing dofs

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/elements/ElementDataEnum.h"
#include "nuto/mechanics/elements/ElementEnum.h"
#include "nuto/mechanics/elements/IpDataEnum.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#else
#include <boost/ptr_container/ptr_vector.hpp>
#endif //ENABLE_SERIALIZATION

#ifdef SHOW_TIME
    #include <ctime>
#endif

#include <iostream>
#include "nuto/base/NuToObject.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/optimize/CallbackHandlerGrid.h"
#include "nuto/optimize/ConjugateGradientGridRed.h"
#include "nuto/mechanics/elements/Brick8N.h"

int main()
{
#define NEIGHBORS
	bool matrixFreeMethod=0; //0 -EBE, 1- NBN, false=0

//	bool matrixFreeMethod=1; //0 -EBE, 1- NBN
	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
	std::cout<<"[NuTo::Grid3D] matrixFreeMethod is";
	if (matrixFreeMethod)
	{
		std::cout<<" NBN \n";
		outputTime<<" NBN  ";

	}
	else
	{
		std::cout<<" EBE \n";
		outputTime<<" EBE  ";
	}

	//double microTomm =0.001;
	double microTomm =1.;

	//   int readFlag = false;
    double PoissonsRatio = 0.2;
    //for local base stiffness matrix
    double YoungsModulus = 1.;

    int rNumVoxel;
    double rVoxelSpacing[3];
    double rGridOrigin[3];
    int rGridDimension[3];


    // create structure
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	// read entries
	NuTo::StructureGrid myGrid(3);
	myGrid.SetVerboseLevel(2);
	// input files have to have a one voxel thick frame without material
	// dimensions are inclusive the frame
	myGrid.ImportFromVtkASCIIFileHeader("InputTest",rGridDimension,rVoxelSpacing,rGridOrigin,rNumVoxel);
//	std::cout<<__FILE__<<"  variab: spac: "<<rVoxelSpacing[0]<< " gridDim: "<<rGridDimension[0]<<std::endl;
	assert (rNumVoxel==rGridDimension[0]*rGridDimension[1]*rGridDimension[2]);
	std::cout<<"[NuTo::Grid3D] numVoxel "<<rNumVoxel<<std::endl;
	outputTime<<rNumVoxel<<"   ";

	std::vector<int> imageValues (rNumVoxel);
	myGrid.ImportFromVtkASCIIFile( "InputTest",imageValues);

	int numGridNodes=(rGridDimension[0]+1)*(rGridDimension[1]+1)*(rGridDimension[2]+1);//all nodes of the grid
    int numBytesEBE=0;

 	//RB
	double Force = -1.;
	bool EnableDisplacementControl = true;
//	bool EnableDisplacementControl = false;
	double BoundaryDisplacement = -1.0;
//	double BoundaryDisplacement = -(rGridDimension[2]*rVoxelSpacing[2]*microTomm)/20.0;
	std::cout<<"[NuTo::Grid3D]  Boundary Displacement: "<<BoundaryDisplacement<<std::endl;
	//calculate one element stiffness matrix with E=1
	NuTo::Structure myHelpStruc(3);

	myHelpStruc.SetVerboseLevel(1);
	// create material law
	int myMat=myHelpStruc.ConstitutiveLawCreate("LinearElastic");
	myHelpStruc.ConstitutiveLawSetPoissonsRatio(myMat, PoissonsRatio);
	myHelpStruc.ConstitutiveLawSetYoungsModulus(myMat, YoungsModulus);

	// create nodes
	NuTo::FullMatrix<double> nodeCoordinates(3, 1);
	NuTo::FullMatrix<int> elementIncidence(8,1);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(0,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(1,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(2,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(3,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(4,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(5,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(6,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm * 0.5;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm * 0.5;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm * 0.5;
	elementIncidence(7,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	// first element create

   // elementIncidence.Info();
	int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
	myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);

	//myHelpStruc.NodeInfo(0);
//	myHelpStruc.NodeBuildGlobalDofs();

	// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	NuTo::FullMatrix<double> stiffnessMatrix;
	NuTo::FullMatrix<int> rows;
	NuTo::FullMatrix<int> coluums;
	myHelpStruc.ElementStiffness(0,stiffnessMatrix,rows,coluums );
//	std::cout<<"Element Stiffness created"<<std::endl;
//	stiffnessMatrix.Info(9,5);
	stiffnessMatrix.WriteToFile("stiffnessMatrix.txt"," ");
	std::vector<double> baseStiffness(24*24);
	for (int i=0;i<24;++i)
	{
		for (int j=0;j<24;++j)
		{
			baseStiffness[(24*i)+j]=stiffnessMatrix(i,j); //row based saved
		}
	}

//	stiffnessMatrix=NuTo::FullMatrix<double>(24,24,baseStiffness);
//	stiffnessMatrix.WriteToFile("baseStiffness.txt"," ");

//	const double rLocalCoordinates[3]={0,0,0};
//	std::vector<double> rDerivativeShapeFunctions(24);
//	NuTo::ElementBase* myElem=myHelpStruc.ElementGetElementPtr(0);
//	NuTo::Brick8N* myBrick=static_cast<NuTo::Brick8N*> (myElem);
//	myBrick->CalculateDerivativeShapeFunctionsLocal(&rLocalCoordinates[0],rDerivativeShapeFunctions);
//	NuTo::FullMatrix<double> derShapeMat(rDerivativeShapeFunctions);
//	derShapeMat.WriteToFile("derivativeShape.txt"," ");


	//grid structure create
	/**************************************************************************/
	//material values are smaller than threshold value
	//color values form 0 255
	//thresholdMaterialValues: 180, 188 ->last value with material
	int thresholdMaterialValue=188; //last value for "air", voxel which does not get element
	//set Modul for each color
	std::vector<double> myMapColorModul(256);
	for(int count=0;count<thresholdMaterialValue;count++)
	{
//		myMapColorModul[count]=1.;
		myMapColorModul[count]=100000.;
	}
	for(int count=thresholdMaterialValue;count<255;count++)
		myMapColorModul[count]=0.;

	size_t numElems=0;
	size_t numNodes=0;
	boost::dynamic_bitset<> nodeExist(numGridNodes); //0 = false, all 0 here
	std::vector<size_t>  edgeId(numNodes);		//saves the id of the edge of this node
	std::vector<int>  nodeId(numGridNodes);		//saves the id of the edges of this node
	std::vector<size_t> voxelId(numElems);		//saves the id of the voxel of this element
	std::vector<double> youngsModulus(0);
	std::vector<int> materialOfElem(numElems);

	// put get routines for a while here
	// get voxel locations of voxels
	std::vector<int> allEdgesAtVoxel(rNumVoxel*8);
//	std::vector<int> allFirstNodesOfElem(numElems);
	for (int element=0;element<rNumVoxel;++element)
	{
//		if (elemExist[element]) // elem exist
//		{
			int numDimxy=element/((rGridDimension[0])*(rGridDimension[1]));
			int numDimx=0;
			int residual1=element%((rGridDimension[0])*(rGridDimension[1]));
			int residual2=0;
			numDimx=residual1/(rGridDimension[0]);
			residual2=residual1%(rGridDimension[0]);
			int rVoxelLocation[3];
			rVoxelLocation[0]=residual2;
			rVoxelLocation[1]=numDimx;
			rVoxelLocation[2]=numDimxy;
			allEdgesAtVoxel[8*element+0] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0];
			allEdgesAtVoxel[8*element+1] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+2] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0] +1;
			allEdgesAtVoxel[8*element+3] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0];

			allEdgesAtVoxel[8*element+4] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0];
			allEdgesAtVoxel[8*element+5] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+6] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allEdgesAtVoxel[8*element+7] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0];

		//	std::cout<<" nodesAtVoxel "<<element<<": "<<allEdgesAtVoxel[8*element+0]<<" "<<allEdgesAtVoxel[8*element+1]<<" "<<allEdgesAtVoxel[8*element+2]<<" "<<allEdgesAtVoxel[8*element+3]<<" "<<allEdgesAtVoxel[8*element+4]<<" "<<allEdgesAtVoxel[8*element+5]<<" "<<allEdgesAtVoxel[8*element+6]<<" "<<allEdgesAtVoxel[8*element+7]<<"\n";
			//also for global vector
			//allVoxelLocation[element]=new int[3];
			//allVoxelLocation[element]=rVoxelLocation;
//		}
//	allFirstNodesOfElem[element]=allEdgesAtVoxel[8*voxelId[element]];

	}

	int numCoeffMat=0;
	bool matExistsAlready=false;
	for (int countVoxels=0;countVoxels<rNumVoxel;++countVoxels)
	{
		if (imageValues[countVoxels]<thresholdMaterialValue)
		{
			voxelId.push_back(countVoxels);
			++numElems;
			for (int node=0;node<8;++node)
			{
				nodeExist.set(allEdgesAtVoxel[8*countVoxels+node],true);
			}

			// nodes are no more set here but not yet elsewhere
			for(int countMat=0;countMat<numCoeffMat;countMat++)
			{
				if (myMapColorModul[imageValues[countVoxels]]==youngsModulus.at(countMat)) //same modulus already used
				{
					materialOfElem.push_back(countMat);
					countMat=numCoeffMat+1;
					matExistsAlready=true;
				}
				else
					matExistsAlready=false;
			}
				if (!matExistsAlready)
				{
					materialOfElem.push_back(numCoeffMat);
					//set youngsModulus and add on material on counter
					youngsModulus.push_back(myMapColorModul[imageValues[countVoxels]]);
					numCoeffMat++;
					matExistsAlready=true; //matrix already added
				}
  		}
	}
    int numDofs=0;
 	for (int node=0;node<numGridNodes;++node)
	{
		if (nodeExist[node])
		{
			edgeId.push_back(node);
		//	std::cout<<" Grid3DRed: edge "<<node<<" node "<<numNodes<<"\n";
			nodeId[node]=numNodes++;
		}
		else
		{
			nodeId[node]=-1;
		}
	}
 	numDofs=3*numNodes;
	std::cout<<"[NuTo::Grid3D] numElems = "<<numElems<<"\n";
	std::cout<<"[NuTo::Grid3D] numNodes = "<<numNodes<<"\n";
	std::cout<<"[NuTo::Grid3D] numDofs = "<<numDofs<<"\n";

	std::cout<<"[NuTo::Grid3D] Young's Modulus = "<<youngsModulus[0]<<"\n";
//	std::cout<<" nodeExist reverse "<<nodeExist<<"\n";
//	std::cout<<" elemExist "<<elemExist<<"\n";

	//NBN**************************************************************************
//	int orderNodeForElem[8][8]={
//		{0,1,4,3,9,10,13,12},
//		{1,2,5,4,10,11,14,13},
//		{4,5,8,7,13,14,17,16},
//		{3,4,7,6,12,13,16,15},
//		{9,10,13,12,18,19,22,21},
//		{10,11,14,13,19,20,23,22},
//		{13,14,17,16,22,23,26,25},
//		{12,13,16,15,21,22,25,24},
//	};
	//field of 27 nodes which contain element and node of element
	std::vector<int> allNodesAtNode;
	std::vector<double> edgeStiffness;
	if(matrixFreeMethod)
	{
		allNodesAtNode.resize(27*numGridNodes);

		int neighborNodes[27];
		neighborNodes[0]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1)-1;
		neighborNodes[1]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1);
		neighborNodes[2]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1)+1;
		neighborNodes[3]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)-1;
		neighborNodes[4]=-(rGridDimension[0]+1)*(rGridDimension[1]+1);
		neighborNodes[5]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)+1;
		neighborNodes[6]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1)-1;
		neighborNodes[7]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1);
		neighborNodes[8]=-(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1)+1;
		neighborNodes[9]=-(rGridDimension[0]+1)-1;
		neighborNodes[10]=-(rGridDimension[0]+1);
		neighborNodes[11]=-(rGridDimension[0]+1)+1;
		neighborNodes[12]=-1;
		neighborNodes[13]=0;
		neighborNodes[14]=+1;
		neighborNodes[15]=+(rGridDimension[0]+1)-1;
		neighborNodes[16]=+(rGridDimension[0]+1);
		neighborNodes[17]=+(rGridDimension[0]+1)+1;
		neighborNodes[18]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1)-1;
		neighborNodes[19]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1);
		neighborNodes[20]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)-(rGridDimension[0]+1)+1;
		neighborNodes[21]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)-1;
		neighborNodes[22]=+(rGridDimension[0]+1)*(rGridDimension[1]+1);
		neighborNodes[23]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)+1;
		neighborNodes[24]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1)-1;
		neighborNodes[25]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1);
		neighborNodes[26]=+(rGridDimension[0]+1)*(rGridDimension[1]+1)+(rGridDimension[0]+1)+1;


	//	std::cout<<"                Neighbors: ";
	//	for (int count=0;count<27;count++)
	//		std::cout<< neighborNodes[count]<<" ";
	//	std::cout<<" "<<std::endl;
	//
		for(int node=0;node<numGridNodes;++node)
		{
			std::vector<int> locNeighbor(27);
			for (int i=0;i<27;++i)
				locNeighbor[i]=node+neighborNodes[i];

	// -------------------------------------------------------------------------------------------
	// special handling of neighbor nodes no longer needed because of the zero margin
	// -------------------------------------------------------------------------------------------
			for(int i=0;i<27;++i)
				allNodesAtNode[27*node+i]=locNeighbor[i];

	//		std::cout<<" Knoten "<<node<< " Neighbors: ";
	//		for (int count=0;count<27;count++)
	//			std::cout<< locNeighbor[count]<<" ";
	//		std::cout<<" "<<std::endl;
		}

		// edgeStiffness 3x3 ************************
		// first field of 8 local nodes of one element
		// second field of 8 local nodes of one element
		// contains NBN ordering local node id
		int orderNode[8][8]={
			{13,14,17,16,22,23,26,25},
			{12,13,16,15,21,22,25,24},
			{9,10,13,12,18,19,22,21},
			{10,11,14,13,19,20,23,22},
			{4,5,8,7,13,14,17,16},
			{3,4,7,6,12,13,16,15},
			{0,1,4,3,9,10,13,12},
			{1,2,5,4,10,11,14,13},
		};

		edgeStiffness.resize(numGridNodes*9*27);

		for (size_t element=0;element<numElems;++element)
		{
			int numMat =materialOfElem[element];
			for (int i=0;i<8;++i)
			{
				for(int j=0;j<8;++j)
				{
					for(int row=0;row<3;++row)
					{
						for(int col=0;col<3;++col)
						{
							edgeStiffness[27*9*allEdgesAtVoxel[i+8*voxelId[element]]+9*orderNode[i][j]+3*row+col]+=youngsModulus[numMat]*baseStiffness[(i*3+row)*24+j*3+col];
						}
					}
				}
			}
		}
//		std::cout<<"\n edgeStiff ";
//		for (int i=0;i<numGridNodes*9*27;++i)
//			std::cout<<edgeStiffness[i]<<" ";
//		std::cout<<"\n";
	}



	//********************************************************************************

	boost::dynamic_bitset<> rDofIsConstraint(3*numNodes); //0 = false, all 0 here

	std::vector<double> displVector(3*numNodes,0.0);// initialized with zero

	//myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");

//---------------------------------------------------------------------------------------//
// boundary conditions
//---------------------------------------------------------------------------------------//
	int numConstraintDofs=0;
// --------------------------------------------------------------------------------//
	//all nodes with x=0
	// set x,y,z direction constraint , zero
//	for (int count = 0;count<numGridNodes;count+= (rGridDimension[0] + 1))
//	{
//		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint x"<< count <<std::endl;
//		rDofIsConstraint.set(count*3+0,true);
//		rDofIsConstraint.set(count*3+1,true);
//		rDofIsConstraint.set(count*3+2,true);
//	}

	// all nodes with y=0 and x=0
	// set y direction zero
//	for (int count =0;count < numGridNodes;count+= (rGridDimension[0] + 1)*(rGridDimension[1] + 1))
//	{
//		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint y "<< count <<std::endl;
//		rDofIsConstraint.set(count*3+1,true);
//	}

	// all nodes with z=0 and x=0
	// set z direction zero
//	for (int count = 0;count<(rGridDimension[0] + 1)*(rGridDimension[1] + 1);count+=(rGridDimension[0] + 1))
//	{
//		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint z "<< count <<std::endl;
//		rDofIsConstraint.set(count*3+2,true);
//	}

	//----------------------------------------------------------------------------------------//
	// Boundary condition,all nodes with x=0 and z=0
	//----------------------------------------------------------------------------------------//
	// set x direction zero
//	for (int count = 0;count<(rGridDimension[0] + 1)*(rGridDimension[1] + 1);count+=(rGridDimension[0] + 1))
//	{
//		rDofIsConstraint.set(count*3,true);
//		++numConstraintDofs;
//	}
	//----------------------------------------------------------------------------------------//
	// Boundary condition,all nodes with x=0
	// Boundary condition,node 0 uy=uz=0 and node = 1*dimx*dimy uz=0
	//----------------------------------------------------------------------------------------//
//	for (int count = 0;count<numGridNodes;count+=(rGridDimension[0] + 1))
//	{
//		rDofIsConstraint.set(count*3,true);
//		++numConstraintDofs;
//	}

//	rDofIsConstraint.set(1,true);
//	++numConstraintDofs;
//	rDofIsConstraint.set(2,true);
//	++numConstraintDofs;
//	rDofIsConstraint.set(3*((rGridDimension[0]+1)*(rGridDimension[1]+1)-1)+2,true);
//	++numConstraintDofs;

	//----------------------------------------------------------------------------------------//

	//----------------------------------------------------------------------------------------//
	// Boundary condition: all nodes with y=0 and z=0
	//Boundary condition:  set y direction zero
////	for (int count = 0;count<(rGridDimension[0] + 1);++count)
//	// Boundary condition: all nodes with y=0
//	std::cout<<"[NuTo::Grid3D] Boundary conditions:  \n";
//	for(int row=0;row<(rGridDimension[2] + 1);++row)
//	{
//		for (int count = 0;count<(rGridDimension[0] + 1);++count)
//		{
//			rDofIsConstraint.set(3*row*(rGridDimension[0] + 1)*(rGridDimension[1] + 1)+count*3+1,true);
//			++numConstraintDofs;
//		}
//
//	}

	//----------------------------------------------------------------------------------------//
	// Boundary condition: all nodes with z=0
	// Boundary condition: set x,y,z direction zero
	//----------------------------------------------------------------------------------------//
	std::cout<<"[NuTo::Grid3D] Boundary conditions: all dofs constraint for z=0 \n";
	for (int count = (rGridDimension[0] + 1)*(rGridDimension[1] + 1);count<2*(rGridDimension[0] + 1)*(rGridDimension[1] + 1);++count)
	{
		if(nodeId[count]>-1)
		{
			rDofIsConstraint.set(nodeId[count]*3+0,true);
			rDofIsConstraint.set(nodeId[count]*3+1,true);
			rDofIsConstraint.set(nodeId[count]*3+2,true);
	//		++numConstraintDofs;
			numConstraintDofs+=3;

		}
	}
	if(numConstraintDofs==0)
		std::cout<<"[NuTo::Grid3D] No boundary conditions set (z=0).  \n";

	//----------------------------------------------------------------------------------------//

	// apply nodes
	std::vector<double>extForces(0);
	if(EnableDisplacementControl)
	{
//		std::cout << "Displacement control" << std::endl;
		// boundary displacments

		//----------------------------------------------------------------------------------------//
		// Boundary condition: all nodes with z=max, uz=BoundaryDisplacement
		std::cout<<"[NuTo::Grid3D] Boundary conditions: Displacement at z=max in z-direction \n";
		int help=numConstraintDofs;
		for (int count = (rGridDimension[0] + 1)*(rGridDimension[1] + 1)*(rGridDimension[2]-1);count<(rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2];++count)
		{
			if(nodeId[count]>-1)
			{
				rDofIsConstraint.set(nodeId[count]*3+2,true);
				displVector[nodeId[count]*3+2]=BoundaryDisplacement;
				++numConstraintDofs;
			}
		}
		if(numConstraintDofs-help==0)
			std::cout<<"[NuTo::Grid3D] No boundary conditions set (z=max).  \n";
		//----------------------------------------------------------------------------------------//

		//----------------------------------------------------------------------------------------//
		// Boundary condition: all nodes with x=max, ux=BoundaryDisplacement
		//----------------------------------------------------------------------------------------//
//		for (int count = rGridDimension[0];count<numGridNodes;count+=(rGridDimension[0] + 1))
//		{
//			rDofIsConstraint.set(count*3,true);
//			displVector[count*3]=BoundaryDisplacement;
//			++numConstraintDofs;
//		}
		//----------------------------------------------------------------------------------------//
	}
	else
	{
		std::cout << "[NuTo::Grid3D] Load control" <<std::endl;
		extForces.resize(3*numGridNodes,0.);
		// Boundary condition: fz=-1
//		for (int count = (rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2];count<numGridNodes;++count)
//		{
//
// 			if(nodeExist[count])
//			{
// 				extForces[count*3+2]=Force;
//			}
//		}
		// for 8 element example, unit force 1
		int count = (rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2];
//		std::cout <<__FILE__<<" "<<__LINE__<< "first node" <<count<<std::endl;
		extForces[count*3+2]=Force/4;
		++count;
		extForces[count*3+2]=Force/2;
		++count;
		extForces[count*3+2]=Force/4;
		++count;
		extForces[count*3+2]=Force/2;
		++count;
		extForces[count*3+2]=Force;
		++count;
		extForces[count*3+2]=Force/2;
		++count;
		extForces[count*3+2]=Force/4;
		++count;
		extForces[count*3+2]=Force/2;
		++count;
		extForces[count*3+2]=Force/4;

		std::cout<<"  extForces ";
		for (int i=0;i<3*numGridNodes;++i)
				std::cout<< extForces[i] <<" ";
			std::cout<<"\n";

	///////////////////////////////////////////////////////////

//		for (int xcount=0;xcount< rGridDimension[0];++xcount)
//		{
//			double nodeForce;
//				if(xCount == 0 || xCount == rGridDimension[0])
//				{
//					nodeForce = Force / (4 *rGridDimension[1] * rGridDimension[2]);
//				}
//				else
//				{
//					nodeForce = Force / (2 *rGridDimension[1] * rGridDimension[2]);
//				}
//
//		}
// 		for (int count = (rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2];count<numGridNodes;++count)
//		{
//
// 			if(nodeExist[count])
//			{
//				if (count==0 || count==(rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2]+rGridDimension[0]
//						||count==(rGridDimension[0] + 1)*(rGridDimension[1] + 1)*rGridDimension[2]+(rGridDimension[0] + 1)*(rGridDimension[1]) || count==(rGridDimension[0] + 1)*(rGridDimension[1] + 1)*(rGridDimension[2]+1)-1 )
//					nodeForce = Force / (4 *rGridDimension[1] * rGridDimension[0]);
//				else if ()
//
//
//				rDofIsConstraint.set(count*3+2,true);
//				displVector[count*3+2]=BoundaryDisplacement;
//				++numConstraintDofs;
//			}
//		}
//
	}
//	std::cout<<"  constraint z ";
//	for (size_t i=0;i<rDofIsConstraint.size();++i)
//			std::cout<< rDofIsConstraint[i] <<" ";
//		std::cout<<"\n";

	outputTime<<numDofs<<"   ";
#ifdef SHOW_TIME
end=clock();
std::cout<<"[NuTo::Grid3D] structure set " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n";
	outputTime<<difftime(end,start)/CLOCKS_PER_SEC<<"   ";
#endif

	std::cout<<"[NuTo::Grid3D] number of dofs "<<numDofs<<" free: "<<numDofs-numConstraintDofs<<" constraint: "<<numConstraintDofs<<"\n";
	// start analysis
//	std::cout<<__FILE__<<" "<<__LINE__<<"  start analysis"<<std::endl;

	NuTo::ConjugateGradientGridRed myOptimizer(numNodes*3);

	myOptimizer.SetVerboseLevel(0);

	myOptimizer.SetParameters(displVector);
	std::cout<<"[NuTo::Grid3D] Parameters set, Anzahl = "<<numDofs<<std::endl;
//	std::cout<<"  voxelId ";
//	for (size_t i=0;i<numElems;++i)
//			std::cout<< voxelId[i] <<" ";
//		std::cout<<"\n";
//
//		if (matrixFreeMethod) //NBN
//	{
//		voxelId.resize(0);
//	}
//	std::cout<<"  nodeId ";
//	for (int i=0;i<numGridNodes;++i)
//			std::cout<< nodeId[i] <<" ";
//		std::cout<<"\n";
//
//	std::cout<<"  edgeId ";
//	for (size_t i=0;i<numNodes;++i)
//			std::cout<< edgeId[i] <<" ";
//		std::cout<<"\n";

	//
	// overhead std constructors
	// with savig neighbors,
	numBytesEBE=7*sizeof(std::vector<double>)+sizeof(rDofIsConstraint);
	numBytesEBE+=sizeof(numNodes)+sizeof(rGridDimension)+sizeof(matrixFreeMethod)
					+sizeof(double)*(
						rDofIsConstraint.num_blocks()+youngsModulus.size()+baseStiffness.size()
						+edgeStiffness.size()+displVector.size()+extForces.size())
				+sizeof(int)*(
						voxelId.size()+edgeId.size()+nodeId.size()+materialOfElem.size()
						+allEdgesAtVoxel.size()+allNodesAtNode.size());
	// with overhead of CGMethode with P: p,r,pr,d,h
	numBytesEBE+=6*sizeof(double)+5*sizeof(int)+5*sizeof(std::vector<double>)+5*sizeof(double)*displVector.size();

	std::cout<<"[NuTo::Grid3D] Required memory: "<<numBytesEBE/(1024.*1024)<<" MiB ,"<<numBytesEBE/(1000000.)<<" MB ("<<numBytesEBE<<" byte) "<<std::endl;
	outputTime<<numBytesEBE/(1000000.)<<"   ";
	outputTime.close();
	myOptimizer.Initialize(numNodes*3,rGridDimension,matrixFreeMethod,voxelId,edgeId,nodeId,
			rDofIsConstraint,youngsModulus,baseStiffness,edgeStiffness,materialOfElem,
			allEdgesAtVoxel,allNodesAtNode,displVector,extForces);
//	myOptimizer.AnsysInput(numGridNodes*3,elemExist,nodeExist,rDofIsConstraint,youngsModulus,rGridDimension,rVoxelSpacing,materialOfElem,allEdgesAtVoxel,displVector);

//	std::cout<<"[NuTo::Grid3D] sizeof std::vector "<<sizeof(std::vector<double>)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof double "<<sizeof(double)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof numNodes "<<sizeof(numNodes)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof rGridDimension"<<" "<<sizeof(rGridDimension)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof matrixFreeMethod "<<sizeof(matrixFreeMethod)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof voxelId"<<" "<<sizeof(voxelId)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof nodeId"<<" "<<sizeof(nodeId)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof rDofIsConstraint "<<8*rDofIsConstraint.num_blocks()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof youndsModulus "<<youngsModulus.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof materialOfElem "<<materialOfElem.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof baseStiffness "<<baseStiffness.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof(edgeStiffness) "<<sizeof(edgeStiffness)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof edgeStiffnes "<<edgeStiffness.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof allNodes "<<allNodesAtNode.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof edgeStiffness.capacity "<<edgeStiffness.capacity()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof nodeExist "<<nodeExist.num_blocks()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof nodeExist "<<nodeExist.bits_per_block<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof nodeExist "<<sizeof(nodeExist)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof nodeExist "<<nodeExist.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof dof "<<rDofIsConstraint.num_blocks()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof dof "<<rDofIsConstraint.bits_per_block<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof dof "<<sizeof(rDofIsConstraint)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof bool "<<sizeof(bool)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof double "<<sizeof(double)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof int "<<sizeof(int)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof displ "<<displVector.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof extForces "<<sizeof(extForces)<<std::endl;
//	std::cout<<"[NuTo::Grid3D] sizeof extForces "<<extForces.size()<<std::endl;
//	std::cout<<"[NuTo::Grid3D] initialized with number of bytes "<<numBytesEBE<<std::endl;
	rDofIsConstraint.clear();
	materialOfElem.clear();
	youngsModulus.clear();
	edgeStiffness.clear();
	allEdgesAtVoxel.clear();
	allNodesAtNode.clear();

	myOptimizer.Optimize();
	myOptimizer.GetParameters(displVector);
	// open file
	std::ofstream file;
//    file.open("displVTK.txt");
//	for(int i=0;i<3*numGridNodes;++i)
//		file<<displVector[i]<<"\n";
//	file.close();

    file.open("displacements.txt");
	for(size_t i=0;i<numNodes;++i)
	{
			file<<displVector[3*i]<<"\n";
			file<<displVector[3*i+1]<<"\n";
			file<<displVector[3*i+2]<<"\n";
	}
	file.close();
	//
	std::vector<double> dispRef;
	std::ifstream input;
	double help=0;
	// result file only with existing nodes
	input.open("result.txt");
	if(input)	// file is open
	{
		int i=0;
		while(!input.eof()) // keep reading untill end-of-file
		{
			input>>help;
			dispRef.push_back(help);
			++i;
		}
		input.close();
		--i; // for last empty line

		if (i==numDofs)
		{
			double squareDiffNorm=0;
			double squareRefNorm=0;
// output of diff and ref only for VTK
//			std::ofstream diffFile;
//			diffFile.open("displDiffVTK.txt");
//			file.open("displRefVTK.txt");
			int k=0;
			for(int i=0;i<numGridNodes;++i)
			{
				if (nodeExist[i])
				{
//					diffFile<<displVector[3*i]-dispRef[3*i]<<"\n";
//					diffFile<<displVector[3*i+1]-dispRef[3*i+1]<<"\n";
//					diffFile<<displVector[3*i+2]-dispRef[3*i+2]<<"\n";
//					file<<dispRef[3*i]<<"\n";
//					file<<dispRef[3*i+1]<<"\n";
//					file<<dispRef[3*i+2]<<"\n";
					squareDiffNorm+=(displVector[3*i]-dispRef[3*k])*(displVector[3*i]-dispRef[3*k]);
					squareDiffNorm+=(displVector[3*i+1]-dispRef[3*k+1])*(displVector[3*i+1]-dispRef[3*k+1]);
					squareDiffNorm+=(displVector[3*i+2]-dispRef[3*k+2])*(displVector[3*i+2]-dispRef[3*k+2]);
					squareRefNorm+=(dispRef[3*k])*(dispRef[3*k]);
					squareRefNorm+=(dispRef[3*k+1])*(dispRef[3*k+1]);
					squareRefNorm+=(dispRef[3*k+2])*(dispRef[3*k+2]);
					++k;
				}
//				else
//				{
//					diffFile<<"0.0 \n";
//					diffFile<<"0.0 \n";
//					diffFile<<"0.0 \n";
//					file<<"0.0 \n";
//					file<<"0.0 \n";
//					file<<"0.0 \n";
//				}
			}
			std::cout<<"[NuTo::Grid3D] squared diff norm " <<squareDiffNorm<<std::endl;
			std::cout<<"[NuTo::Grid3D] error " <<sqrt(squareDiffNorm)/sqrt(squareRefNorm)*100<<" %"<<std::endl;
		}
		else
			std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (wrong size).\n";

	}
	else
		std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (no result file).\n";
	return 0;
}
