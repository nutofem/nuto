// $Id$
//! @author Andrea Keßler, ISM
//! @date September 2012
//! @brief ... standard class for multigrid routines

#include "nuto/optimize/MultiGrid.h"

void NuTo::MultiGrid::SetStructure(NuTo::StructureGrid* rpStructureHandler)
{
	mpStructureHandler=rpStructureHandler;
}

void NuTo::MultiGrid::SetParameters(std::vector<double>& u )
{
	mpStructureHandler->SetParameters(u);
}
std::vector<double>&  NuTo::MultiGrid::GetParameters()
{
	return (mpStructureHandler->GetParameters());
}

int NuTo::MultiGrid::Initialize()
{
	std::string restrictionType="TWENTYSEVENPOINTS";
	std::cout<<"[MultiGrid] ("<<__LINE__<<") Initialize "<<mpStructureHandler<<"\n";
    assert(mpStructureHandler->GetDimension()==3);
    std::cout<<"[MultiGrid] "<<mpStructureHandler->GetDimension()<<"\n";

	mpStructureHandler->CalculateMultiGridCorrolations(restrictionType,mRestriction,mProlongation);
	// for 27 points prolongation is equal restriction * 1/8.
	//ProlongationMatrix();

    assert(mpStructureHandler->GetGridDimension()[2]);


    std::vector<size_t> rGridDimension= mpStructureHandler->GetGridDimension();
    std::cout<<"[MultiGrid] (line "<<__LINE__<<") dim "<<rGridDimension[0]<<"\n";
    assert(mpStructureHandler->GetDimension()==3);
    int dim=mpStructureHandler->GetDimension();

	double mem=sizeof(mpStructureHandler);


	while(rGridDimension[0]%2==0 && rGridDimension[1]%2==0 && rGridDimension[2]%2==0  &&
		  rGridDimension[0]>4 && rGridDimension[1]>4 && rGridDimension[2]>4) // delete this part later
	{
		NuTo::StructureGrid* rCoarseGrid=new NuTo::StructureGrid(dim);
		if(mpGridPtr.empty())
		{
			mpStructureHandler->SetCoarserGridLevel(rCoarseGrid);
		}
		else
		{
			std::cout<<"[MultiGrid] coarse grid number "<<mNumGrids<<" create.\n";
			NuTo::StructureGrid* rGrid=GetGridPtr(mNumGrids-2);
			rGrid->SetCoarserGridLevel(rCoarseGrid);
		}
		std::cout<<"[MultiGrid] coarse grid number "<<mNumGrids<<" created.\n";
		++mNumGrids;
		mpGridPtr.push_back(rCoarseGrid);
		mem+=sizeof(rCoarseGrid);
		std::cout<<"[MultiGrid] mpGridPtr "<<mpGridPtr.size()<<" .\n";

		rGridDimension=rCoarseGrid->GetGridDimension();
	}

	// Cycle type ///
	std::string cycleType="VCYCLE";
	mCycle.resize(2*(mNumGrids+mNumGrids-2));

	for(int i=0;i<mNumGrids-1;++i)
	{
		mCycle[2*i]=i;
		mCycle[2*i+1]=1;
	}

	mCycle[2*(mNumGrids-1)]=mNumGrids-1;
	mCycle[2*(mNumGrids-1)+1]=-1;
	int j=2;
	for(int i=mNumGrids;i<(int) (mCycle.size()/2);++i)
	{
		mCycle[2*i]=mNumGrids-j++;
		mCycle[2*i+1]=-1;
	}
	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
	outputTime<<mem<<"  ";
	outputTime.close();
	std::cout<<"[MultiGrid3D] memory "<<mem/1000.<<" MB \n";

	return 0;
}


int NuTo::MultiGrid::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startMG,endMG;
    startMG=clock();
#endif
	std::cout<<"[MultiGrid] Optimize \n";
	// start with coarse grid as start solution
//	SetCurrentGridNumber(mNumGrids);

	//name of future routine
	double minObjective=1e-6;	//square of value taken
	int cycle=0;
	double residual=1.;
	const int rMaxGradientCalls=5;

	std::cout<<"[MultiGrid] "<<__LINE__<<" mNumGrids "<<mNumGrids<<"\n";

	// new //
	assert(mCycle.size()%2==0);
	std::cout<<" cycle ";
	for(size_t i=0;i<mCycle.size();++i)
	{
		std::cout<<mCycle[i]<<" ";
	}
	std::cout<<"\n";
	std::cout<<"mCycle.size() "<<mCycle.size()<<"\n";
	bool button=true;

	while (button)
	{
		NuTo::StructureGrid* rGrid;
		for(size_t step=0;step<mCycle.size()/2;++step)
		{
			if(mCycle[2*step]==0)
				rGrid=mpStructureHandler;
			else
				rGrid=GetGridPtr(mCycle[2*step]-1);

			if(mVerboseLevel>2)

			{
				std::cout<<"[MultiGrid] ("<<__LINE__<<") optimize grid "<<mCycle[2*step]<<"\n";
				std::cout<<"[MultiGrid] "<<" u ";
				for(int i=0;i<rGrid->GetNumNodes()*3;++i)
				{
					std::cout<<rGrid->GetParameters()[i]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" r ";
				for(int i=0;i<rGrid->GetNumNodes()*3;++i)
				{
					std::cout<<rGrid->GetResidual()[i]<<" ";
				}
				std::cout<<"\n";
			}
		//	SolveOneVCycle();
			// set jacobi for all grids instead of the coarsest which has to be solved not smoothed
			if (mCycle[2*step]+1!=mNumGrids)
			{
				NuTo::Jacobi rOptiJacobi(rGrid->GetNumNodes()*3);
				rOptiJacobi.SetCallback(rGrid);
				rOptiJacobi.SetVerboseLevel(0);
//				if(mCycle[2*step+1]==-1)
//				if(cycle==0 && mCycle[2*step+1]==1)
//					rOptiJacobi.SetMaxGradientCalls(1);
//				else
					rOptiJacobi.SetMaxGradientCalls((int) rGrid->GetNumNodes()*0.01);
//					rOptiJacobi.SetMaxGradientCalls(10);
//					rOptimizer.SetParameters(rGrid->GetParameters());
//				else
				rOptiJacobi.Optimize();
				std::cout<<"[MultiGrid] smoothed \n";
				if(mCycle[2*step]==0 && cycle!=0)
				{
					residual=rOptiJacobi.GetObjective();
					std::cout<<"[MultiGrid] ("<<__LINE__<<") residual "<<residual<<"\n";
				}
			}
			else // solve coarsest grid with CG
			{
				NuTo::ConjugateGradientGrid rOptiCG(rGrid->GetNumNodes()*3);
				rOptiCG.SetCallback(rGrid);
				rOptiCG.SetVerboseLevel(0);
				rOptiCG.Optimize();
			}
			// pre-smoothing
			std::cout<<"[MultiGrid] ("<<__LINE__<<") Grid "<<mCycle[2*step]<<" optimized \n";
			// set break criteria for while-clauses
			if(residual<minObjective)
			{
				button=false;
				step=mCycle.size()/2;
			}

			if(mVerboseLevel>2)
			{
				std::cout<<"[MultiGrid] "<<" u ";
				for(int i=0;i<rGrid->GetNumNodes()*3;++i)
				{
					std::cout<<rGrid->GetParameters()[i]<<" ";
				};
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" r ";
				for(int i=0;i<rGrid->GetNumNodes()*3;++i)
				{
					std::cout<<rGrid->GetResidual()[i]<<" ";
				}
				std::cout<<"\n";
			}
			 //Restriction value equal 1
			if(mCycle[2*step+1]==1) //Restriction
			{
//				if(cycle==0)
//					rGrid->StartRestriction(mRestriction);
//				else
					rGrid->Restriction(mRestriction);
			}
			else if (mCycle[2*step+1]==-1)
				rGrid->Prolongation(mProlongation); // on coarse grid -> put paras to fine grid
//#ifdef SHOW_TIME
//			if(mCycle[2*step]==0)
//			{
//				if (mShowTime)
//					std::cout<<"[MultiGrid] " << difftime(clock(),startMG)/CLOCKS_PER_SEC << "sec" << std::endl;
//			}
//#endif //SHOW_TIME
		}
		++cycle;
		std::cout<<"[MultiGrid] ("<<__LINE__<<") Residual "<<residual<<" Cycle "<<cycle<<" \n";
	}
	std::cout<<"[MultiGrid] Optimized in "<<cycle<<" Cycles. \n";

#ifdef SHOW_TIME
    endMG=clock();
    if (mShowTime)
        std::cout<<"[MultiGrid] " << difftime(endMG,startMG)/CLOCKS_PER_SEC << "sec" << std::endl;
	std::fstream outputTime;
	std::string filename = "timeOutput";
    outputTime.open(filename,std::fstream::out|std::fstream::app);
 	outputTime<<(difftime(endMG,startMG)/CLOCKS_PER_SEC)<<"   "<<cycle<<"\n";
	outputTime.close();

#endif

	return 0;
}

void NuTo::MultiGrid::Info()const
{
	throw OptimizeException(std::string("[MultiGrid] not implemented."));
}

//! @brief gives number of grids
//! @return number of grids
int NuTo::MultiGrid::GetNumGrids()
{
	return mNumGrids;
}
//! @brief get pointer to  grid
//! @param grid identifier
//! @return pointer to grid
NuTo::StructureGrid* NuTo::MultiGrid::GetGridPtr(int rIdent)
{
	return &mpGridPtr.at(rIdent);
}

//! @brief set current grid number
//! @param rGridNumber ... grid number
void NuTo::MultiGrid::SetCurrentGridNumber(int rGridNumber)
{
	mCurrentGridNumber=rGridNumber;
}

//! @brief gives number of relaxation iterations - smoothing on fine grid
//! @return number of relaxation iterations
int NuTo::MultiGrid::GetNumRelaxIterations()
{
	return mNumRelaxIterations;
}

//! @brief set number of relaxation iterations - smoothing on fine grid
//! @param numIterations ... number of relaxation iterations
void NuTo::MultiGrid::SetNumRelaxIterations(int numIterations)
{
	mNumRelaxIterations=numIterations;
}

