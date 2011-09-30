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
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/mechanics/elements/Brick8N.h"

int main()
{
	//double microTomm =0.001;
	double microTomm =1;

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
	myGrid.ImportFromVtkASCIIFileHeader("InputTest",rGridDimension,rVoxelSpacing,rGridOrigin,rNumVoxel);
	std::cout<<__FILE__<<"  variab: spac: "<<rVoxelSpacing[0]<< " gridDim: "<<rGridDimension[0]<<std::endl;
	assert (rNumVoxel==rGridDimension[0]*rGridDimension[1]*rGridDimension[2]);
	std::cout<<__FILE__<<" "<<__LINE__<<" numVoxel "<<rNumVoxel<<std::endl;

	std::vector<int> imageValues (rNumVoxel);
	myGrid.ImportFromVtkASCIIFile( "InputTest",imageValues);
//	std::cout<<"image values ";
//	for(int i=0;i<rNumVoxel;++i)
//		std::cout<<imageValues[i] <<" ";
//	std::cout<<"\n";

	int numGridNodes=(rGridDimension[0]+1)*(rGridDimension[1]+1)*(rGridDimension[2]+1);//all nodes of the grid

	//RB
	//double Force = 1.;
	bool EnableDisplacementControl = true;
	double BoundaryDisplacement = 1.0;
//	double BoundaryDisplacement = (rGridDimension[0]*rVoxelSpacing[0]*microTomm)/20.0;
	std::cout<<__FILE__<<"  Boundary Displacement: "<<BoundaryDisplacement<<std::endl;
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

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(0,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(1,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(2,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = -rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(3,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(4,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = -rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(5,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(6,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	nodeCoordinates(0, 0) = -rVoxelSpacing[0] * microTomm / 2;
	nodeCoordinates(1, 0) = rVoxelSpacing[1] * microTomm / 2;
	nodeCoordinates(2, 0) = rVoxelSpacing[2] * microTomm / 2;
	elementIncidence(7,0)=myHelpStruc.NodeCreate("displacements", nodeCoordinates);

	// first element create

   // elementIncidence.Info();
	int myHelpElement=myHelpStruc.ElementCreate("Brick8N", elementIncidence);
	//myHelpStruc.ElementSetConstitutiveLaw(element,"Material1");
	myHelpStruc.ElementSetConstitutiveLaw(myHelpElement,myMat);

	//myHelpStruc.NodeInfo(0);
	myHelpStruc.NodeBuildGlobalDofs();

	// build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
	NuTo::FullMatrix<double> stiffnessMatrix;
	NuTo::FullMatrix<int> rows;
	NuTo::FullMatrix<int> coluums;
	myHelpStruc.ElementStiffness(0,stiffnessMatrix,rows,coluums );
	std::vector<double> baseStiffness(24*24);
	for (int i=0;i<24;++i)
	{
		for (int j=0;j<24;++j)
		{
// Test variant
//			baseStiffness[(24*i)+j]=1.; //row based saved
			baseStiffness[(24*i)+j]=stiffnessMatrix(i,j); //row based saved
		}
	}
	std::cout<<"Element Stiffness created"<<std::endl;
//	stiffnessMatrix.Info(12,6);
	stiffnessMatrix.WriteToFile("stiffnessMatrix.txt"," ");
//	stiffnessMatrix=baseStiffness;
//	stiffnessMatrix.WriteToFile("baseStiffness.txt"," ");

	const double rLocalCoordinates[3]={0,0,0};
	std::vector<double> rDerivativeShapeFunctions(24);
	NuTo::ElementBase* myElem=myHelpStruc.ElementGetElementPtr(0);
	NuTo::Brick8N* myBrick=static_cast<NuTo::Brick8N*> (myElem);
	myBrick->CalculateDerivativeShapeFunctionsLocal(&rLocalCoordinates[0],rDerivativeShapeFunctions);

	NuTo::FullMatrix<double> derShapeMat(rDerivativeShapeFunctions);
	derShapeMat.WriteToFile("derivativeShape.txt"," ");


	//grid structure create
	/**************************************************************************/
	//material values are smaller than threshold value
	//color values form 0 255
	//thresholdMaterialValues: 180, 188 ->last value with material
	int thresholdMaterialValue=188; //last value for "air", voxel which does not get element
	//set Modul for each color
	std::vector<double> myMapColorModul(256);
	for(int count=0;count<thresholdMaterialValue;count++)
		myMapColorModul[count]=100000.;
	for(int count=thresholdMaterialValue;count<255;count++)
		myMapColorModul[count]=0.;

	std::cout<<" Young's Modulus = 100000\n";
	boost::dynamic_bitset<> nodeExist(numGridNodes); //0 = false, all 0 here
	boost::dynamic_bitset<> elemExist(rNumVoxel); //0 = false, all 0 here
	std::vector<double> youngsModulus(0);
	std::vector<int> materialOfElem(rNumVoxel);

	//std::vector<int*> allVoxelLocation(0,rNumVoxel);
//		std::vector<int> allFirstNodes(-1,rNumVoxel);
	std::vector<int> allNodesAtVoxel(rNumVoxel*8);


	// put get routines for a while here
	// get voxel locations of voxels
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
			allNodesAtVoxel[8*element+0] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0];
			allNodesAtVoxel[8*element+1] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allNodesAtVoxel[8*element+2] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0] +1;
			allNodesAtVoxel[8*element+3] = rVoxelLocation[2]*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0];

			allNodesAtVoxel[8*element+4] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0];
			allNodesAtVoxel[8*element+5] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + rVoxelLocation[1]     * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allNodesAtVoxel[8*element+6] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0]+1;
			allNodesAtVoxel[8*element+7] = (rVoxelLocation[2]+1)*(rGridDimension[0]+1)*(rGridDimension[1]+1) + (rVoxelLocation[1]+1) * (rGridDimension[1]+1) + rVoxelLocation[0];

//			std::cout<<" nodesAtVoxel "<<element<<": "<<allNodesAtVoxel[8*element+0]<<" "<<allNodesAtVoxel[8*element+1]<<" "<<allNodesAtVoxel[8*element+2]<<" "<<allNodesAtVoxel[8*element+3]<<" "<<allNodesAtVoxel[8*element+4]<<" "<<allNodesAtVoxel[8*element+5]<<" "<<allNodesAtVoxel[8*element+6]<<" "<<allNodesAtVoxel[8*element+7]<<"\n";
			//allFirstNodes[element]=allNodesAtVoxel[8*element];
			//also for global vector
			//allVoxelLocation[element]=new int[3];
			//allVoxelLocation[element]=rVoxelLocation;
//		}
	}
	//myGrid.SetAllElementIds();
	//myGrid.SetAllNodeIdsAtNode();
	//myGrid.SetAllPartCoefficientMatrix0();
	int numCoeffMat=0;
	bool matExistsAlready=false;
	for (int countVoxels=0;countVoxels<rNumVoxel;++countVoxels)
	{
		if (imageValues[countVoxels]<thresholdMaterialValue)
		{
			elemExist.set(countVoxels,true);
			for (int node=0;node<8;++node)
			{
				nodeExist.set(allNodesAtVoxel[8*countVoxels+node],true);
			}
			for(int countMat=0;countMat<numCoeffMat;countMat++)
			{
				if (myMapColorModul[imageValues[countVoxels]]==youngsModulus.at(countMat)) //same modulus already used
				{
					materialOfElem[countVoxels]=countMat;
					countMat=numCoeffMat+1;
					matExistsAlready=true;
				}
				else
					matExistsAlready=false;
			}
				if (!matExistsAlready)
				{
					materialOfElem[countVoxels]=numCoeffMat;
					//set youngsModulus and add on material on counter
					youngsModulus.push_back(myMapColorModul[imageValues[countVoxels]]);
					numCoeffMat++;
					matExistsAlready=true; //matrix already added
				}
  		}
 		else
 		{
// 			std::cout<<"elem number "<<countVoxels<<" does not exist \n";
 			materialOfElem[countVoxels]=-1;
 		}

	}
//	std::cout<<"Grid created"<<std::endl;
//	std::cout<<" nodeExist "<<nodeExist<<"\n";
	// set dofs
	boost::dynamic_bitset<> rDofIsConstraint(3*numGridNodes); //0 = false, all 0 here
//	std::cout<<"  constraint orig ";
//	for (int i=0;i<3*numGridNodes;++i)
//			std::cout<< rDofIsConstraint[i] <<" ";
//		std::cout<<"\n";
	for (int i=0;i<numGridNodes;++i)
	{
		if (!nodeExist[i]) //node does not exist
		{
			// dof is constraint
			rDofIsConstraint.flip(3*i);
			rDofIsConstraint.flip(3*i+1);
			rDofIsConstraint.flip(3*i+2);
		}
	}

	std::vector<double> displVector(3*numGridNodes,0.0);// initialized with zero

	//myMapColorModul.WriteToFile("$HOME/develop/nuto/MapColorModul.txt"," ");
	//generiert Knoten, Freiheitsgrade noch nicht gesetzt

	// boundary conditions
	int NumElementsX = rGridDimension[0];
	int NumElementsY = rGridDimension[1];
	int NumElementsZ = rGridDimension[2];

	NuTo::FullMatrix<double> direction(3,1);
	direction(0,0)= 1;
	direction(1,0)= 0;
	direction(2,0)= 0;

	for (int count = 0;count<numGridNodes;count+= (NumElementsX + 1))
	{
		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint x"<< count <<std::endl;
		rDofIsConstraint.set(count*3+0,true);
	}

	direction(0,0)= 0;
	direction(1,0)= 1;
	direction(2,0)= 0;
	// Test Randbedingungen y=0
	// for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);++count)
	for (int count =0;count < numGridNodes;count+= (NumElementsX + 1)*(NumElementsY + 1))
	{
		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint y "<< count <<std::endl;
		rDofIsConstraint.set(count*3+1,true);
	}

	direction(0,0)= 0;
	direction(1,0)= 0;
	direction(2,0)= 1;
	// Test Randbedingungen z=0
	//for (int count =0;count < (NumElementsX + 1)*(NumElementsY + 1)*(NumElementsZ + 1);++count)
	for (int count = 0;count<(NumElementsX + 1)*(NumElementsY + 1);count+=(NumElementsX + 1))
	{
		//std::cout<<__FILE__<<" "<<__LINE__<<" node constraint z "<< count <<std::endl;
		rDofIsConstraint.set(count*3+2,true);
	}

	// apply nodes
	if(EnableDisplacementControl)
	{
//		std::cout << "Displacement control" << std::endl;
		// boundary displacments
		direction(0,0)= 1;
		direction(1,0)= 0;
		direction(2,0)= 0;
		NuTo::FullMatrix<double> displacements(3,1);
		displacements(0,0)= BoundaryDisplacement;
		displacements(1,0)= 0;
		displacements(2,0)= 0;

		for(int zCount = 0; zCount < NumElementsZ + 1; ++zCount)
		{
//			std::cout << zCount << std::endl;
			for(int yCount = 0; yCount < NumElementsY + 1; ++yCount)
			{
				int node = zCount * (NumElementsX + 1) * (NumElementsY + 1) + yCount * (NumElementsX + 1) + NumElementsX;
//				 std::cout << "node displacement constraint: " << node << std::endl;
				if (nodeExist[node])
				{
					displVector[node*3]=BoundaryDisplacement;
					rDofIsConstraint.set(node*3,true);
				}
			}
		}
	}
	else
	{
		std::cout <<__FILE__<<" "<<__LINE__<< "Load control" << "not implemented"<<std::endl;
		return 1;
	}
//	std::cout<<"  constraint ";
//	for (int i=0;i<3*numGridNodes;++i)
//			std::cout<< rDofIsConstraint[i] <<" ";
//		std::cout<<"\n";

#ifdef SHOW_TIME
end=clock();
std::cout<<"[NuTo::StructureGrid3D] structure set " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif

	// start analysis
	std::cout<<__FILE__<<" "<<__LINE__<<"  start analysis"<<std::endl;
	// build global dof numbering
	myGrid.SetVerboseLevel(0);

	//NuTo::ConjugateGradientGrid myOptimizer((unsigned int) numGridNodes*3);
	NuTo::ConjugateGradientGrid myOptimizer((unsigned int) numGridNodes*3);
	std::cout<<__FILE__<<" "<<__LINE__<<"  optimizer created"<<std::endl;

	std::cout<<__FILE__<<" "<<__LINE__<<"  "<<std::endl;
	myOptimizer.SetVerboseLevel(1);

	myOptimizer.SetParameters(displVector);
	std::cout<<__FILE__<<" "<<__LINE__<<"  Parameters set, Anzahl "<<myOptimizer.GetNumParameters()<<std::endl;
//#ifdef ENABLE_MECHANICS
//		myOptimizer.SetGridStructure(&myGrid);
// #endif // ENABLE_MECHANICS

	//set callback routines for the calculation of the objective function, gradient etc
	//this works, because Neural network has been derived from CallbackHandler of the optimization module
	myOptimizer.Initialize(numGridNodes*3,elemExist,nodeExist,rDofIsConstraint,youngsModulus,baseStiffness,materialOfElem,allNodesAtVoxel,displVector);
	myOptimizer.AnsysInput(numGridNodes*3,elemExist,nodeExist,rDofIsConstraint,youngsModulus,rGridDimension,rVoxelSpacing,materialOfElem,allNodesAtVoxel,displVector);
	std::cout<<__FILE__<<" "<<__LINE__<<"  initialized"<<std::endl;
	myOptimizer.Optimize();
	myOptimizer.GetParameters(displVector);
	NuTo::FullMatrix<double> displacements(displVector);
	displacements.WriteToFile("displacements.txt"," ");
	//
	NuTo::FullMatrix<double> dispAnsys(numGridNodes*3,1);
	dispAnsys.ReadFromFile("/ismhome/tests/ansys/result.txt");
	if(numGridNodes*3==dispAnsys.GetNumColumns())
	{
		displacements-=dispAnsys;
		displacements.WriteToFile("displDiff.txt"," ");
		std::cout<<" diff norm " <<displacements.Norm()<<std::endl;
		std::cout<<" error " <<displacements.Norm()/dispAnsys.Norm()*100<<" %"<<std::endl;
	}
	else
		std::cout<<" Comparision with Ansys results is not possible.\n";
	return 0;
}
