// $Id$
//! @author Andrea Keßler, ISM
//! @date September 2012
//! @brief ... standard class for multigrid routines

#include "nuto/optimize/MultiGrid.h"

void NuTo::MultiGrid::SetStructure(NuTo::StructureGrid* rpStructureHandler)
{
	mpStructureHandler=rpStructureHandler;
}

int NuTo::MultiGrid::Initialize()
{
	std::string restrictionType="TWENTYSEVENPOINTS";
	mpStructureHandler->CalculateMultiGridCorrolations(restrictionType,mRestriction,mProlongation);
	// for 27 points prolonagation is equal restriction * 1/8.
	//ProlongationMatrix();

    std::vector<size_t> rGridDimension= mpStructureHandler->GetGridDimension();
    std::cout<<__LINE__<<" "<<"dim "<<rGridDimension[0]<<"\n";
    int dim=mpStructureHandler->GetDimension();
	while(rGridDimension[0]%2==0 && rGridDimension[1]%2==0 && rGridDimension[2]%2==0  &&
		  rGridDimension[0]>4 && rGridDimension[1]>4 && rGridDimension[2]>4) // delete this part later


	{
		NuTo::StructureGrid* rCoarseGrid=new NuTo::StructureGrid(dim);
		if(mpGridPtr.empty())
			mpStructureHandler->SetCoarserGridLevel(rCoarseGrid);
		else
		{
			std::cout<<"[MultiGrid] coarse grid number "<<mNumGrids<<" create.\n";
			NuTo::StructureGrid* rGrid=GetGridPtr(mNumGrids-2);
			rGrid->SetCoarserGridLevel(rCoarseGrid);
		}
		std::cout<<"[MultiGrid] coarse grid number "<<mNumGrids<<" created.\n";
		++mNumGrids;
		mpGridPtr.push_back(rCoarseGrid);
		std::cout<<"[MultiGrid] mpGridPtr "<<mpGridPtr.size()<<" .\n";


		rGridDimension=rCoarseGrid->GetGridDimension();
	}
	return 0;
}


int NuTo::MultiGrid::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startMG,endMG;
    startMG=clock();
#endif

	// start with coarse grid as start solution
//	SetCurrentGridNumber(mNumGrids);

	//name of future routine
	double minObjective=1e-12;	//square of value taken
	int cycle=0;
	double residual=1.;
	int i=0,j=1;
	const int rMaxGradientCalls=2;

	std::cout<<"[MultiGrid] "<<__LINE__<<" mNumGrids "<<mNumGrids<<"\n";
	int precision = 7;

	while (residual> minObjective)
	{
		NuTo::StructureGrid* rGrid;
		while(i<mNumGrids)
		{
			if(i==0)
			{
				if(cycle==0)
					rGrid=mpStructureHandler;
				else
				{
					++i;++j;
					rGrid=GetGridPtr(i-1);
				}

			}
			else
				rGrid=GetGridPtr(i-1);

		//	SolveOneVCycle();
			NuTo::ConjugateGradientGrid rOptimizer(rGrid->GetNumNodes()*3);
			rOptimizer.SetCallback(rGrid);
			rOptimizer.SetMaxGradientCalls(rMaxGradientCalls);
			std::vector<double> rParameters(rGrid->GetNumNodes()*3);
			if(mVerboseLevel>2)
			{
				std::cout<<" [MultiGrid] optimize grid "<<i<<"\n";
//			std::cout<<"[MultiGrid] "<<" u ";
//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
//			{
//				std::cout<<rGrid->GetParameters()[i]<<" ";
//			}
//			std::cout<<"\n";
//			std::cout<<"[MultiGrid] "<<" r ";
//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
//			{
//				std::cout<<rGrid->GetResidual()[i]<<" ";
//			}
//			std::cout<<"\n";
				std::cout.precision(precision);
						//std::setw(width);
				std::cout<<"[MultiGrid] "<<" ux ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" uy ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i+1]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" uz ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i+2]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" rx ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" ry ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i+1]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" rz ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i+2]<<" ";
					}
					std::cout<<"\n";
			}
			rOptimizer.Optimize();
				std::cout<<" [MultiGrid] Grid "<<i<<" optimized \n";
			if(mVerboseLevel>2)
			{
	//			std::cout<<"[MultiGrid] "<<" u ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetParameters()[i]<<" ";
	//			};
	//			std::cout<<"\n";
	//			std::cout<<"[MultiGrid] "<<" r ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetResidual()[i]<<" ";
	//			}
	//			std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" ux ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" uy ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i+1]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" uz ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetParameters()[3*i+2]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" rx ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" ry ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i+1]<<" ";
					}
					std::cout<<"\n";
					std::cout<<"[MultiGrid] "<<" rz ";
					for(int i=0;i<rGrid->GetNumNodes();++i)
					{
						std::cout<<rGrid->GetResidual()[3*i+2]<<" ";
					}
					std::cout<<"\n";
			}
			if(i<j)
				rGrid->Restriction(mRestriction);
			else
				rGrid->Prolongation(mProlongation); // on coarse grid -> put paras to fine grid
			++i;
			if(j<mNumGrids-1)
				++j;
			else
				j=i-2;
		}
		--i;
		while(i>0)
		{
			--i;
			if(j>0)
				--j;
			else
				j=i+1;
			NuTo::StructureGrid* rGrid=NULL;
			if(i==0)
				rGrid=mpStructureHandler;
			else
				rGrid=GetGridPtr(i-1);
			NuTo::ConjugateGradientGrid rOptimizer(rGrid->GetNumNodes()*3);
			rOptimizer.SetCallback(rGrid);
			rOptimizer.SetMaxGradientCalls(rMaxGradientCalls);
			std::vector<double> rParameters(rGrid->GetNumNodes()*3);
			if(mVerboseLevel>2)
			{
				std::cout<<" [MultiGrid] optimize grid "<<i<<"\n";
	//			std::cout<<"[MultiGrid] "<<" u ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetParameters()[i]<<" ";
	//			};
	//			std::cout<<"\n";
	//			std::cout<<"[MultiGrid] "<<" r ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetResidual()[i]<<" ";
	//			}
	//			std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" ux ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" uy ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i+1]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" uz ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i+2]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" rx ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" ry ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i+1]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" rz ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i+2]<<" ";
				}
				std::cout<<"\n";
			}
			rOptimizer.Optimize();
				std::cout<<" [MultiGrid] Grid "<<i<<" optimized \n";
			if(mVerboseLevel>2)
			{
	//			std::cout<<"[MultiGrid] "<<" u ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetParameters()[i]<<" ";
	//			};
	//			std::cout<<"\n";
	//			std::cout<<"[MultiGrid] "<<" r ";
	//			for(int i=0;i<rGrid->GetNumNodes()*3;++i)
	//			{
	//				std::cout<<rGrid->GetResidual()[i]<<" ";
	//			}
	//			std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" ux ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" uy ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i+1]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" uz ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetParameters()[3*i+2]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" rx ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" ry ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i+1]<<" ";
				}
				std::cout<<"\n";
				std::cout<<"[MultiGrid] "<<" rz ";
				for(int i=0;i<rGrid->GetNumNodes();++i)
				{
					std::cout<<rGrid->GetResidual()[3*i+2]<<" ";
				}
				std::cout<<"\n";
			}
			residual=rOptimizer.GetObjective();

			if(i>j)
				rGrid->Prolongation(mProlongation); // on coarse grid -> put paras to fine grid
			else
				rGrid->Restriction(mRestriction);
		}
		// end one cycle
		++cycle;
//		std::cout<<"[MultiGrid] residual "<<sqrt(residual)<<" \n";
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

