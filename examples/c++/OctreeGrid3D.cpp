// $Id$
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
#include "nuto/mechanics/structures/grid/OctreeGrid.h"
#include "nuto/mechanics/structures/grid/MultiGridStructure.h"
#ifdef ENABLE_OPTIMIZE
#include "nuto/optimize/ConjugateGradientGrid.h"
#include "nuto/optimize/MisesWielandt.h"
#endif

int main()
{
    bool matrixFreeMethod=0; //0 -EBE, 1- NBN, false=0
//	bool matrixFreeMethod=1; //0 -EBE, 1- NBN


//    double PoissonsRatio = 0.;
    double PoissonsRatio = 0.2;
    // create structure
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	// read entries
	NuTo::OctreeGrid myGrid(3); // also creates CallbackHandler
	myGrid.SetVerboseLevel(0);
//	std::string inputFile="InputTest";
	std::string inputFile="InputStructureGrid3D";
	myGrid.ImportFromVtkASCIIFileHeader(inputFile);

	//calculate one element stiffness matrix with E=1

	myGrid.SetMatrixFreeMethod(matrixFreeMethod);
	myGrid.SetBasisElementStiffnessMatrix(PoissonsRatio,0);

	//grid structure create
	/**************************************************************************/
	//material values are smaller than threshold value
	//color values form 0 255
	//thresholdMaterialValues: 180, 188 ->last value with material
	int thresholdMaterialValue=188; //last value for "air", voxel which does not get element
	std::vector<double> myMapColorModul(256);

	//set Modul for each color
	double MaterialYoungsModulus=100000.;
	for(int count=0;count<thresholdMaterialValue;count++)
//		myMapColorModul[count]=1.;
		myMapColorModul[count]=MaterialYoungsModulus;
	for(int count=thresholdMaterialValue;count<255;count++)
		myMapColorModul[count]=0.;

	// create grid structure
//	myGrid.SetMaxOctreeLevels(1);
	myGrid.CreateOctree(thresholdMaterialValue,inputFile,myMapColorModul);
	size_t numNodes=myGrid.GetNumNodes();

	//----------------------------------------------------------------------------------------//
	// Boundary condition: x=0, ux=0; y=0, uy=0; z=0,uz=0;z0max,uz=-1
	//----------------------------------------------------------------------------------------//

	//	bool EnableDisplacementControl = false;
	double BoundaryDisplacement = 1.0;
//	double BoundaryDisplacement = (double)(myGrid.GetGridDimension()[1])*0.01;
	// diplacement vector plus one dof for calculation with non existing neighbor dofs
	std::vector<double> rDisplVector(3*(3*numNodes+1),0.0);// initialized with zero
	//for z=0,x,y -all ux=0
	const std::vector<size_t> rGridDimension=myGrid.GetGridDimension();
	// Attention: no considering of frame nodes for octree
	double rValue=0.;
	// rGridLocation ... region of constrained nodes xmin,xmax,ymin,ymax,zmin,zmax
	size_t rGridLocation[6]={0};

	//x=0, ux=0
	size_t direction=0;
	rGridLocation[0]=0;
	rGridLocation[1]=0;
	rGridLocation[2]=0;
	rGridLocation[3]=rGridDimension[1]+1;
	rGridLocation[4]=0;
	rGridLocation[5]=rGridDimension[2]+1;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

	//y=0; uy=0
	direction=1;
	rGridLocation[0]=0;
	rGridLocation[1]=rGridDimension[0]+1;
	rGridLocation[2]=0;
	rGridLocation[3]=0;
	rGridLocation[4]=0;
	rGridLocation[5]=rGridDimension[2]+1;

	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,-1,rDisplVector);


	//  uz=0; for one plain
	direction=2;
	rGridLocation[0]=0;
	rGridLocation[1]=rGridDimension[0]+1;
	rGridLocation[2]=0;
	rGridLocation[3]=rGridDimension[1]+1;
	rGridLocation[4]=0;
	rGridLocation[5]=0;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,rValue,rDisplVector);

//y=max, uy=-1
	direction=1;
	rGridLocation[0]=0;
	rGridLocation[1]=rGridDimension[0]+1;
	rGridLocation[2]=rGridDimension[1];
	rGridLocation[3]=rGridDimension[1]+1;
	rGridLocation[4]=0;
	rGridLocation[5]=rGridDimension[2]+1;
	rValue=BoundaryDisplacement;
	myGrid.SetDisplacementConstraints(direction,rGridLocation,1,rDisplVector);
//	myGrid.SetDisplacementConstraints(direction,rGridLocation,0,rDisplVector);

	std::cout<<"rDisplVector  ";
	for(size_t i=0;i<numNodes*3;++i)
		std::cout<<rDisplVector[i]<<" ";
	std::cout<<"\n";

	myGrid.ExportVTKUnstructuredGridDataFile("./outputFileGeo.vtk");
	myGrid.AnsysInput(rDisplVector);
#ifdef SHOW_TIME
end=clock();
#endif
	double mem=sizeof(myGrid);
	myGrid.Info();
	std::cout<<"-   Poissons ratio .............................. "<<PoissonsRatio<<"\n\n";
	std::cout<<"Boundary displacement .......................... "<<BoundaryDisplacement<<"\n";
	std::cout<<"Direction of external force/displacement ....... z\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";
	std::cout<<"Time for structure initialization............... "<<difftime(end,start)/CLOCKS_PER_SEC<<"(sec)\n";
	std::cout<<"Allocated memory ............................... "<<mem/1000.<<"(MB)\n";
	std::cout<<"-----------------------------------------------------------------------------------\n";

#ifdef ENABLE_OPTIMIZE
	enum SolMethod // enum for solution method
	{
		MG, 	// multigrid method
		JCG, 	// jacobi conjugate gradient method
		MGCG, 	//multigrid preconditioned conjugate gradient method
	} solMeth=JCG;

	if(solMeth==MG)
	{

//		NuTo::MultiGridStructure myMultiGrid;
//		myMultiGrid.SetVerboseLevel(0);
//		myMultiGrid.SetStructure(&myGrid);
//		myMultiGrid.SetUseMultiGridAsPreconditoner(false);
//		myMultiGrid.SetMaxCycle(numDofs/10);
//		myMultiGrid.SetNumPreSmoothingSteps(1);
//		myMultiGrid.SetNumPostSmoothingSteps(1);
//
//		myMultiGrid.Initialize();
//		myMultiGrid.Info();
	}
	else if(solMeth==JCG)
	{
		// start analysis
		NuTo::ConjugateGradientGrid myOptimizer(numNodes*3);
		myOptimizer.SetVerboseLevel(5);
		myGrid.SetMisesWielandt(false);
	//	myGrid.SetWeightingFactor(1);

		myOptimizer.SetCallback( (&myGrid));
		myOptimizer.Info();
		myOptimizer.Optimize();
	}

#else //ENABLE_OPTIMIZE
	std::cout<<"[OctreeGrid3D] Solution is not possible. Module optimize is not loaded.\n";
#endif //ENABLE_OPTIMIZE
	rDisplVector=myGrid.GetParameters();
	//one last hanging node correction on u
	myGrid.HangingNodesCorrection(rDisplVector);

	std::ofstream file;
//	int precision = 15;
//	std::cout.precision(precision);
    file.open("displacements.txt");
	for(size_t i=0;i<numNodes;++i)
	{
			file<<rDisplVector[3*i]<<"\n";
			file<<rDisplVector[3*i+1]<<"\n";
			file<<rDisplVector[3*i+2]<<"\n";
	}
	file.close();
	std::vector<double> dispRef;
	std::ifstream input;
	double help=0;
	// result file only with existing nodes
//	input.open("result.txt");
//	if(input)	// file is open
//	{
//		size_t count=0;
//		while(!input.eof()) // keep reading untill end-of-file
//		{
//			input>>help;
//			dispRef.push_back(help);
//			++count;
//		}
//		input.close();
//		--count; // for last empty line
//
//		if (count==numNodes*3)
//		{
//			double squareDiffNorm=0;
//			double squareRefNorm=0;
//// output of diff and ref only for VTK
////			std::ofstream diffFile;
////			diffFile.open("displDiffVTK.txt");
////			file.open("displRefVTK.txt");
//			for(size_t i=0;i<numNodes;++i)
//			{
////					diffFile<<displVector[3*i]-dispRef[3*i]<<"\n";
////					diffFile<<displVector[3*i+1]-dispRef[3*i+1]<<"\n";
////					diffFile<<displVector[3*i+2]-dispRef[3*i+2]<<"\n";
////					file<<dispRef[3*i]<<"\n";
////					file<<dispRef[3*i+1]<<"\n";
////					file<<dispRef[3*i+2]<<"\n";
//				squareDiffNorm+=(rDisplVector[3*i]-dispRef[3*i])*(rDisplVector[3*i]-dispRef[3*i]);
//				squareDiffNorm+=(rDisplVector[3*i+1]-dispRef[3*i+1])*(rDisplVector[3*i+1]-dispRef[3*i+1]);
//				squareDiffNorm+=(rDisplVector[3*i+2]-dispRef[3*i+2])*(rDisplVector[3*i+2]-dispRef[3*i+2]);
//				squareRefNorm+=(dispRef[3*i])*(dispRef[3*i]);
//				squareRefNorm+=(dispRef[3*i+1])*(dispRef[3*i+1]);
//				squareRefNorm+=(dispRef[3*i+2])*(dispRef[3*i+2]);
//			}
//		std::cout<<"[NuTo::Grid3D] squared diff norm " <<squareDiffNorm<<std::endl;
//		std::cout<<"[NuTo::Grid3D] error " <<sqrt(squareDiffNorm)/sqrt(squareRefNorm)*100<<" %"<<std::endl;
//		}
//		else
//			std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (wrong size).\n";
//	}
//	else
//		std::cout<<"[NuTo::Grid3D] Comparison with reference results is not possible (no result file).\n";

	myGrid.ExportVTKUnstructuredGridDataFile("./OctreeResults.vtk");


	std::cout<<"OctreeGrid3D END\n";
	return 0;
}
