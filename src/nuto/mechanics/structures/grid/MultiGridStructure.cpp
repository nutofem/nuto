// $Id $
//! @author Andrea Keßler, ISM
//! @date March 2013
//! @brief ... standard class for MultiGridStructure structure

#include "nuto/mechanics/structures/grid/MultiGridStructure.h"

//#include "nuto/optimize/Jacobi.h"
//#include "nuto/optimize/ConjugateGradientGrid.h"

void NuTo::MultiGridStructure::SetStructure(NuTo::StructureGrid* rpStructureHandler)
{
	mpStructureHandler=rpStructureHandler;
//	if(mpGridPtr.empty())
//		mpGridPtr.push_back(rpStructureHandler);
}

void NuTo::MultiGridStructure::Hessian(std::vector<double>&  rDiagHessian)
{
	if(mUseAsPreconditioner)
	{
		std::vector<double> rParameters(rDiagHessian.size(),0.);
		MultiGridSolve(rParameters,rDiagHessian);
		for(size_t i=0;i<rDiagHessian.size();++i)
			rDiagHessian[i]=rParameters[i];
	}
	else
		mpStructureHandler->HessianDiag(rDiagHessian);
}
void NuTo::MultiGridStructure::Gradient (std::vector<double>& rValues,std::vector<double>& rGradient)
{
	if(!mpStructureHandler->GetMatrixFreeMethod())
		mpStructureHandler->CalculateMatrixVectorProductEBE(rValues,rGradient);
	else
		mpStructureHandler->CalculateMatrixVectorProductNBN(rValues,rGradient);
}

std::vector<double>&  NuTo::MultiGridStructure::GetParameters()
{
	return mpStructureHandler->GetParameters();
}
std::vector<double>&  NuTo::MultiGridStructure::GetRightHandSide()
{
	return mpStructureHandler->GetRightHandSide();
}

void NuTo::MultiGridStructure::SetParameters(std::vector<double>& rParameters)
{
	mpStructureHandler->SetParameters(rParameters);
}

void NuTo::MultiGridStructure::SetRightHandSide(std::vector<double>& rRightHandSide)
{
	mpStructureHandler->SetRightHandSide(rRightHandSide);
}

int NuTo::MultiGridStructure::Initialize()
{
	std::string restrictionType="TWENTYSEVENPOINTS";
	//! @TODO change back or unifrom
   assert(mpStructureHandler->GetDimension()==3);

	mpStructureHandler->CalculateMultiGridCorrolations(restrictionType,mRestriction,mProlongation);
	// for 27 points prolongation is equal restriction * 1/8.
	//ProlongationMatrix();

    assert(mpStructureHandler->GetGridDimension()[2]);


    std::vector<size_t> rGridDimension= mpStructureHandler->GetGridDimension();

    assert(mpStructureHandler->GetDimension()==3);
    int dim=mpStructureHandler->GetDimension();

	double mem=sizeof(mpStructureHandler);

	while(rGridDimension[0]%2==0 && rGridDimension[1]%2==0 && rGridDimension[2]%2==0
			&&rGridDimension[0]>4 && rGridDimension[1]>4 && rGridDimension[2]>4 // delete this part later
			)
	{
		NuTo::StructureGrid* rCoarseGrid=new NuTo::StructureGrid(dim);
		NuTo::StructureGrid* rGrid=mpStructureHandler;
		if(!mpGridPtr.empty())
			rGrid=GetGridPtr(mNumGrids-2);
		rGrid->SetCoarserGridLevel(rCoarseGrid);
		++mNumGrids;
		mpGridPtr.push_back(rCoarseGrid);
		mem+=sizeof(rCoarseGrid);

		rGridDimension=rCoarseGrid->GetGridDimension();
	}

	// Cycle type ///
	std::string cycleType="VCYCLE";
//	std::string cycleType="WCYCLE";
//	std::string cycleType="VWCYCLE";
//	mCycle.resize(2*(mNumGrids+mNumGrids-2));
//
//	for(int i=0;i<mNumGrids-1;++i)
//	{
//		mCycle[2*i]=i;
//		mCycle[2*i+1]=1;
//	}
//
//	mCycle[2*(mNumGrids-1)]=mNumGrids-1;
//	mCycle[2*(mNumGrids-1)+1]=-1;
//	int j=2;
//	for(int i=mNumGrids;i<(int) (mCycle.size()/2);++i)
//	{
//		mCycle[2*i]=mNumGrids-j++;
//		mCycle[2*i+1]=-1;
//	}
	mCycle.resize(mNumGrids);
	if(cycleType=="VCYCLE")
	{
		for(int i=0;i<mNumGrids;++i)
			mCycle[i]=1;
	}
	else if(cycleType=="WCYCLE")
	{
		mCycle[0]=1;
		mCycle[mNumGrids-1]=1;
		for(int i=1;i<mNumGrids-1;++i)
			mCycle[i]=2;
	}
	else if(cycleType=="FMVCYCLE")
	{
		mCycle[0]=0;
		for(int i=1;i<mNumGrids;++i)
			mCycle[i]=1;
	}
	return 0;
}
int NuTo::MultiGridStructure::MultiGridSolve(std::vector<double>& rSolution,
		std::vector<double>& rRightHandSide)
{
#ifdef SHOW_TIME
    std::clock_t startMG,endMG;
    startMG=clock();
#endif
	int cycle=0;
	double kappa=1.;
//	std::cout<<"MultiGridSolve\n";
	NuTo::StructureGrid* rGrid=mpStructureHandler;
	size_t numPara=rGrid->GetNumNodes()*3;
	double rAccuracySquare=1e-6*1e-6;
//	double rAccuracySquare=1e-6/(double) numPara;
	double rIntialNorm=0.;
	std::vector<double> rResidual(rSolution.size(),0.);
	std::vector<double> rError(rSolution.size(),0.);


	rGrid->Gradient(rSolution,rResidual);

	for(size_t i=0;i<numPara;++i)
	{
		rResidual[i]*=-1;
		rResidual[i]+=rRightHandSide[i];
		rIntialNorm+=(rResidual[i]*rResidual[i]);
	}
	// --------------------------------------------------------------------------
//	std::cout<<"\n start residual "<<__LINE__<<"\n";
//	std::cout <<" mg r ";
//	for(size_t i=1;i<numPara;i+=3)
//		std::cout<<rResidual[i]<<" ";
//	std::cout<<"\n";
//	// --------------------------------------------------------------------------
//	std::cout<<" rIntialNorm "<<rIntialNorm<<"\n";


	while(++cycle<mMaxCycles && kappa>rAccuracySquare)
	{
		rError.assign((numPara+1)*3,0.);
		if(mCycle[0]==0)//FMV- cycle
		{
			std::cout<<"Attention: add restiction to the level.\n";
			for(int rGridLevel=mNumGrids-1;rGridLevel>0;--rGridLevel)
				GridLevelSolve(rGridLevel,rError,rResidual);
		}
		else
			GridLevelSolve(0,rError,rResidual);
//			GridLevelSolve(0,rSolution,rRightHandSide);

//		// update parameters
		for(size_t i=0;i<numPara;++i)
			rSolution[i]+=rError[i];

		rGrid->Gradient(rSolution,rResidual);
		for(size_t i=0;i<numPara;++i)
		{
			rResidual[i]*=-1;
			rResidual[i]+=rRightHandSide[i];
		}
//		// --------------------------------------------------------------------------
//		std::cout<<"\n next rhs "<<__LINE__<<"\n";
//		std::cout <<" mg r ";
//		for(size_t i=1;i<numPara;i+=3)
//			std::cout<<rResidual[i]<<" ";
//		std::cout<<"\n";
//		std::cout <<" mg v ";
//		for(size_t i=1;i<numPara;i+=3)
//			std::cout<<rSolution[i]<<" ";
//		std::cout<<"\n";
//		// --------------------------------------------------------------------------

		kappa=0.;
		for(size_t i=0;i<numPara;++i)
			kappa+=(rResidual[i]*rResidual[i]);

		// relative norm
//		rAccuracySquare*=rIntialNorm;
		kappa/=rIntialNorm;
//		std::cout<<"MultiGridSolve relative error "<<kappa<<" vs. norm "<<rAccuracySquare<<"\n";
	}
#ifdef SHOW_TIME
    endMG=clock();
    if (mShowTime && mVerboseLevel>0)
		std::cout<< "[MultiGridSolve] Elapsed time (sec)............. " << difftime(endMG,startMG)/CLOCKS_PER_SEC << std::endl;
#endif

	int returnValue=0;
	if (kappa<rAccuracySquare)
		returnValue = NORMGRADIENT;
	else if(cycle>=mMaxCycles)
		returnValue = MAXCYCLES;

	if (mVerboseLevel>0)
	{
		std::cout<<"MultiGridSolve relative error "<<kappa<<" vs. norm "<<rAccuracySquare<<"\n";
		std::cout<< " "  << std::endl;
		std::cout<< "[MultiGridSolve] "  << std::endl;
		std::cout<< "Number of Cylces............. " << cycle << std::endl;
		std::cout<< "Active convergence criterion..... " ;
		switch (returnValue)
		{
			case MAXCYCLES:
				std::cout<< "Maximum number of cycles reached." << std::endl;
				break;
			case NORMGRADIENT:
				std::cout<< "Norm of residual smaller than prescribed value." << std::endl;
//				std::cout<< "Relative Norm of residual to intial residual smaller than prescribed value." << std::endl;
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;
		int width=10;
		int precision=6;
		std::cout.precision(precision);
		std::cout << std::setw(width)<< "[MultiGridSolve] displacements " ;

		for (size_t count=0; count<numPara; count++)
		{
			std::cout << std::setw(width)<<rSolution[count] << "   " ;

		}
		std::cout << std::endl;

	}

	return returnValue;
}

void NuTo::MultiGridStructure::GridLevelSolve(int rGridLevel,std::vector<double>& rParameters,
		std::vector<double>& rRightHandSide)
{
#ifdef ENABLE_OPTIMIZE
	int cylcetype=0;
	NuTo::StructureGrid* rGrid;
	if(rGridLevel==0)
		rGrid=mpStructureHandler;
	else
		rGrid=GetGridPtr(rGridLevel-1);

	std::vector<double> rResidual(rParameters.size(),0.);

	// to avoid infinite loop
	rGrid->SetUseDiagHessian(true);

	while(mCycle[rGridLevel]>cylcetype++)
	{
		size_t numPara=rGrid->GetNumNodes()*3;

		if(rGridLevel==mNumGrids-1) // coarsest grid - exact solution
		{
			NuTo::ConjugateGradientGrid rOptiCG(numPara);
			rOptiCG.SetCallback(rGrid);
			// CG takes parameter and rightHandSide from grid structure
			rGrid->SetRightHandSide(rRightHandSide);
			rGrid->SetParameters(rParameters);
			rOptiCG.SetVerboseLevel(0);
			rOptiCG.SetMaxGradientCalls(numPara);
			rOptiCG.Optimize();
			// needed?
			rParameters=rGrid->GetParameters();
		}
		else
		{
			NuTo::Jacobi rOptiJacobi(numPara);
			rOptiJacobi.SetCallback(rGrid);
			rOptiJacobi.SetVerboseLevel(0);
			rOptiJacobi.SetMaxGradientCalls(mNumPreSmoothingSteps);
			rOptiJacobi.Optimize(rParameters,rRightHandSide);
			// get residual of smoothed error or parameters
			// compute residual =rhs - Kv
			rGrid->Gradient(rParameters,rResidual);

			//res=rhs of smoothed e
			for(size_t i=0;i<numPara;++i)
			{
				rResidual[i]*=-1.0;
				rResidual[i]+=rRightHandSide[i];
			}

			//restriction
			std::vector<double> rResidualCoarser;// set size in restriction
			// compute corase rhs from residual
			rGrid->Restriction(mRestriction,rResidual,rResidualCoarser);

			std::vector<double> rParametersCoarser(rResidualCoarser.size(),0.);
			//corase grid correction/solution Ke=r
			GridLevelSolve(rGridLevel+1,rParametersCoarser,rResidualCoarser);
			//prolongaton and correction
			// v <- v +v coarse
			rGrid->Prolongation(mProlongation,rParametersCoarser,rParameters); // on coarse grid -> put paras to fine grid
			// --------------------------------------------------------------------------
//			std::cout<<"\n mg level "<<rGridLevel<<" before smoothing; Line: "<<__LINE__<<"\n";
//			std::cout<<" mg r ";
//			for(size_t i=1;i<numPara;i+=3)
//				std::cout<<rResidual[i]<<" ";
//			std::cout<<"\n";
//			std::cout<<" mg u ";
//			for(size_t i=1;i<numPara;i+=3)
//				std::cout<<rParameters[i]<<" ";
//			std::cout<<"\n";
			// --------------------------------------------------------------------------

			rOptiJacobi.SetMaxGradientCalls(mNumPostSmoothingSteps);
			rOptiJacobi.Optimize(rParameters,rRightHandSide);
		}
	}
	// to reset initial boolean
	rGrid->SetUseDiagHessian(false);
#else //ENABLE_OPTIMIZE
	std::cout<<"\n GridlevelSolve needs enabled Module Optimize.\n";
#endif //ENABLE_OPTIMIZE
}
void NuTo::MultiGridStructure::Info()const
{
	std::cout<<"\nMultiGrid Info\n "
			"-----------------------------------------------------------------------------------\n";
	std::cout<<"Number of grids ................................ " <<mNumGrids<<"\n";
	std::cout<<"Restriction type ............................... " <<mRestrictionType<<"\n";
	std::cout<<"Number of pre smoothing steps .................. " <<mNumPreSmoothingSteps<<"\n";
	std::cout<<"Number of post smoothing steps ................. " <<mNumPostSmoothingSteps<<"\n";
	std::cout<<"Use as preconditioner .......................... " <<mUseAsPreconditioner<<"\n";
	std::cout<<"Multigrid cycle 1 -VCycle, 2-WCycle ............ ";
	for(size_t i=0;i<mCycle.size();++i)		std::cout<<mCycle[i]<<" ";
	std::cout<<"\n";
	std::cout<<"Maximal number of cycles ....................... " <<mMaxCycles<<"\n";
	std::cout<<"-----------------------------------------------------------------------------------\n\n";
}

//! @brief gives number of grids
//! @return number of grids
int NuTo::MultiGridStructure::GetNumGrids()
{
	return mNumGrids;
}
//! @brief get pointer to  grid
//! @param grid identifier
//! @return pointer to grid
NuTo::StructureGrid* NuTo::MultiGridStructure::GetGridPtr(int rIdent)
{
	return &mpGridPtr.at(rIdent);
}

//! @brief set current grid number
//! @param rGridNumber ... grid number
void NuTo::MultiGridStructure::SetCurrentGridNumber(int rGridNumber)
{
	mCurrentGridNumber=rGridNumber;
}

//! @brief gives number of relaxation iterations - smoothing on fine grid
 //! @return number of relaxation iterations
 int  NuTo::MultiGridStructure::GetNumPreSmoothingSteps()
 {
 	return mNumPreSmoothingSteps;
 }
//! @brief set number of relaxation iterations - smoothing on fine grid
//! @param numIterations ... number of relaxation iterations
void NuTo::MultiGridStructure::SetNumPreSmoothingSteps(int rSteps)
{
	mNumPreSmoothingSteps=rSteps;
}

//! @brief gives number of relaxation iterations - smoothing on fine grid
 //! @return number of relaxation iterations
 int  NuTo::MultiGridStructure::GetNumPostSmoothingSteps()
 {
 	return mNumPostSmoothingSteps;
 }
//! @brief set number of relaxation iterations - smoothing on fine grid
//! @param numIterations ... number of relaxation iterations
void NuTo::MultiGridStructure::SetNumPostSmoothingSteps(int rSteps)
{
	mNumPostSmoothingSteps=rSteps;
}

void NuTo::MultiGridStructure::SetUseMultiGridAsPreconditoner(bool rUseAsPreconditioner)
{
	mUseAsPreconditioner=rUseAsPreconditioner;
}
bool NuTo::MultiGridStructure::GetUseMultiGridAsPreconditoner()
{
	return mUseAsPreconditioner;
}
void NuTo::MultiGridStructure::SetMaxCycle(int rMaxCycles)
{
	mMaxCycles=rMaxCycles;
}
int NuTo::MultiGridStructure::GetMaxCycle()
{
	return mMaxCycles;
}
// export to Vtk Datafile
void NuTo::MultiGridStructure::ExportVTKStructuredDataFile(int rGridLevel, const std::string& rFilename)
{
	std::ofstream file(rFilename.c_str());
	if (!file.is_open())
		throw MechanicsException(std::string("[NuTo::MultiGridStructure::ExportVTKStructuredDataFile] Error opening file ")+rFilename.c_str());

	// header /////////////////////////////////////////////////////////////////
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "Data file was generated by NuTo" << std::endl;
	file << "ASCII" << std::endl;
	file << "DATASET STRUCTURED_POINTS" << std::endl;

	size_t countNodes=0;

	StructureGrid* rGrid=NULL;
	if(rGridLevel==0)
	{
		rGrid=mpStructureHandler;
		size_t numNodes=rGrid->GetNumNodes();
		size_t numEdges=(rGrid->mGridDimension[0]+1)*(rGrid->mGridDimension[1]+1)*(rGrid->mGridDimension[2]+1);
		///////////////////////////////////////////////////////////////////////////
		file << "DIMENSIONS "<<(rGrid->mGridDimension[0])+1<<" "<<(rGrid->mGridDimension[1])+1<<" "<<(rGrid->mGridDimension[2])+1<<"\n";
		file << "SPACING "<<rGrid->GetVoxelSpacing()[0]<<" "<<rGrid->GetVoxelSpacing()[1]<<" "<<rGrid->GetVoxelSpacing()[2]<<"\n";
		file << "ORIGIN "<<rGrid->GetGridOrigin()[0]<<" "<<rGrid->GetGridOrigin()[1]<<" "<<rGrid->GetGridOrigin()[2]<<"\n";
		file << "CELL_DATA "<<rGrid->GetNumVoxels()<<"\n";
		file << "SCALARS imageData int 1\n";
		file << "LOOKUP_TABLE default\n";
		using namespace boost::spirit::classic;
	   std::ifstream input(rGrid->mImageDataFile, std::ios::in);
		if (input.is_open() == false)
			throw MechanicsException("[MultiGridStructure::ExportVTKStructuredDataFile] error opening input file.");
		std::string line;
		for (int count=0;count<10;count++)
			getline (input, line);

		int value;
		size_t count=0;
		size_t numVoxels=rGrid->GetNumVoxels();
		while(getline(input,line))
		{
			std::istringstream iss(line);
			while(iss >> value && count<numVoxels)
			{
				file << value<<" \n" ;
				++count;
			}
		}
		 // close file
		input.close();
		file << "POINT_DATA "<<numEdges<<"\n";
		file << "VECTORS displacements double \n";
	//			   std::cout<<" VTK mDisplacements "<<mDisplacements.size()<<"\n";
		countNodes=0;
		for (size_t i=0;i<numEdges;++i) // loop over grid points
		{
			if(rGrid->mNodeId[i]<numNodes)
			{
				file<<rGrid->GetParameters()[3*rGrid->mNodeId[i]]<<" " <<rGrid->GetParameters()[3*rGrid->mNodeId[i]+1]<<" "<<rGrid->GetParameters()[3*rGrid->mNodeId[i]+2]<<"\n";
				++countNodes;
			}
			else
				file<<"0.0 0.0 0.0\n";
		}
		assert(countNodes==numNodes);

	}
	else
	{
		rGrid=GetGridPtr(rGridLevel-1);

		file << "DIMENSIONS "<<rGrid->mGridDimension[0]+1<<" "<<rGrid->mGridDimension[1]+1<<" "<<rGrid->mGridDimension[2]+1<<"\n";
		file << "SPACING "<<rGrid->GetVoxelSpacing()[0]<<" "<<rGrid->GetVoxelSpacing()[1]<<" "<<rGrid->GetVoxelSpacing()[2]<<"\n";
		file << "ORIGIN "<<rGrid->GetGridOrigin()[0]<<" "<<rGrid->GetGridOrigin()[1]<<" "<<rGrid->GetGridOrigin()[2]<<"\n";
//		file << "CELL_DATA "<<rGrid->GetNumVoxels()<<"\n";
//		file << "SCALARS coarser_grid_"<<rGridLevel<<" int 1\n";
//		file << "LOOKUP_TABLE default\n";
		size_t numVoxels=rGrid->GetNumVoxels();
//		size_t numNodes=rGrid->GetNumNodes();
//		size_t numEdges=(rGrid->mGridDimension[0]+1)*(rGrid->mGridDimension[1]+1)*(rGrid->mGridDimension[2]+1);

		size_t element=0;
//		for(size_t voxel=0;voxel<numVoxels;++voxel)
//		{
//			if((rGrid->mVoxelId[element])==voxel)
//			{
//				// only visualization with one material so far
//				file<<"0\n";
//				++element;
//			}
//			else
//				file<<"255\n";
//		}
		file << "CELL_DATA "<<rGrid->GetNumVoxels()<<"\n";
		file << "SCALARS YoungsModulus"<<rGridLevel<<" double\n";
		file << "LOOKUP_TABLE default\n";
		element=0;
		for(size_t voxel=0;voxel<numVoxels;++voxel)
		{
			if((rGrid->mVoxelId[element])==voxel)
			{
				// only visualization with one material so far
				file<<rGrid->mYoungsModulus[element]<<"\n";
				++element;
			}
			else
				file<<"0.\n";
		}
	}
//	size_t numNodes=rGrid->GetNumNodes();
//	size_t numEdges=(rGrid->mGridDimension[0]+1)*(rGrid->mGridDimension[1]+1)*(rGrid->mGridDimension[2]+1);
//
//	countNodes=0;
//	file << "POINT_DATA "<<numEdges<<"\n";
//	file << "VECTORS constraints int \n";
//	for (size_t i=0;i<numEdges;++i) // loop over grid points
//	{
//		if(rGrid->mNodeId[i]<numNodes)
//		{
//			file<<rGrid->mDofIsConstraint[3*rGrid->mNodeId[i]]<<" " <<rGrid->mDofIsConstraint[3*rGrid->mNodeId[i]+1]<<" "<<rGrid->mDofIsConstraint[3*rGrid->mNodeId[i]+2]<<"\n";
//			++countNodes;
//		}
//		else
//			file<<"5 5 5\n";
//	}
//	assert(countNodes==numNodes);
    file.close();
}
