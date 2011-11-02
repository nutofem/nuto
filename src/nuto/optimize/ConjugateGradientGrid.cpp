// $Id $

#include "nuto/optimize/ConjugateGradientGrid.h"
#define machine_precision 1e-15
//sqrt machine_precision

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
#define tol 1e-8

int NuTo::ConjugateGradientGrid::Optimize()
{
#ifdef SHOW_TIME
    std::clock_t startOpt,endOpt;
    startOpt=clock();
#endif
	double alpha=0,
		   beta=0,
		   alphaNumerator=0,
		   alphaDenominator=0,
		   betaNumerator=0,
		   squaredNorm=0;

    int numFunctionCalls(0),   // number of function calls
		 numGradientCalls(0),   // number of gradient calls
		 numHessianCalls(0),    // number of hessian calls
		 curIteration(0),       //number of iterations
		 curCycle(0);           //number of iterations without restart

 	if (mVerboseLevel>5)
	{
		std::cout<<__FILE__<<" "<<__LINE__<<" mMatrixFreeMethod is";
		if (mMatrixFreeMethod)
		{
			std::cout<<" NBN \n";
		}
		else
			std::cout<<" EBE \n";
	}

	optimization_return_attributes returnValue;

	// u  = parameters
	// r  = gradient
	// pr = preconditioned gradient = p r
	// d  = search direction
	// h  = scaled search direction = K d

#ifdef ENABLE_EIGEN
	std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] Eigen vectors are used.\n";
	myType u;
	u.MapAligned(mParameters.data(),mNumParameters);

	if(extForces.size()==0)
		myType f(mNumParameters);
	else
		myType f(mForces.data(),mNumParameters);

	myType r(mNumParameters);
	myType pr;
	myType p(mNumParameters);
	myType h(mNumParameters);
	myType d(mNumParameters);
#else
	//std::unique_ptr
//	myType u(&parameters[0]);
//
//	myType r(new double [mNumParameters]);
//	myType pr;
//	myType p(new double [mNumParameters]);
//	myType h(new double [mNumParameters]);
//	myType d(new double [mNumParameters]);
//	myType u(&parameters[0]);

	//std::vector
	std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] std::vector is used.\n";
	myType &u=mParameters;
	myType &f=mForces;
	myType r(mNumParameters);
	myType pr(mNumParameters);
	myType p(mNumParameters,1);
	myType h(mNumParameters);
	myType d(mNumParameters);
#endif
	int precision = 6;
	int width = 10;

	if (mVerboseLevel>2)
		std::cout<< __FILE__<<" "<<__LINE__<< " Para "<< mNumParameters << std::endl;
	bool converged(false);
	double mAccuracyGradientScaled = mAccuracyGradient;
	std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] gradient accuracy "<<mAccuracyGradientScaled <<std::endl;

	int localMaxGradientCalls=mNumParameters;
//	int localMaxGradientCalls=2*mNumParameters;
	if (localMaxGradientCalls<mMaxGradientCalls)
		SetMaxGradientCalls(localMaxGradientCalls);

	// calculate objective
	numFunctionCalls++;
	if (numFunctionCalls>mMaxFunctionCalls)
	{
		converged = true;
		returnValue = MAXHESSIANCALLS;
	}

	// calculate Diag preonditioner;
//	HessianDiag(p);

//	reduce iterations with accepting an error - for later
//	CalcScalingFactors(numHessianCalls,p);

//	 print p
//	std::cout<<" p ";
//	for(int i=0;i<mNumParameters;i+=3)
//		std::cout<<p[i]<<" ";
//	std::cout<<"\n";


	while(!converged)
	{
		numGradientCalls++;
		 if (numGradientCalls>mMaxGradientCalls)
		 {
			 converged = true;
			 returnValue = MAXGRADIENTCALLS;
			 break;
			 //return MAXGRADIENTCALLS;
		 }

//determine gradient

		// all  mNumParameters repeat start equations
		//if (curCycle%mNumParameters==0)
		// no upate with start equations
		if (curCycle==0)
		{
//			CalcScalingFactors(numHessianCalls,p);
			if (numHessianCalls>mMaxHessianCalls)
			{
				converged = true;
				returnValue = MAXHESSIANCALLS;
				break;
			}
			//calculate gradient as a start solution
			if (mVerboseLevel>3)
				std::cout<<__FILE__<<" "<<__LINE__<<" calc start direction"<<std::endl;
			if(!mMatrixFreeMethod)
				CalculateMatrixVectorProductEBE(u,r);
			else
				CalculateMatrixVectorProductNBN(u,r);

			//get global external force vector
			if (mVerboseLevel>3)
			{
				if(mForces.empty())
				std::cout<<" no external force.\n";
			else
				std::cout<<" with external force.\n";
			}

#ifdef ENABLE_EIGEN
			r*=-1;
			if(!mForces.empty())
				r+=f;
			//or r=KU-f test
//			r-=f;
#else
			for(int i=0;i<mNumParameters;++i)
			{
				if(r[i]!=0)
					r[i]*=-1;
				if(!mForces.empty())
					r[i]+=f[i];
			}
#endif
			// print r
//			std::cout<<" r ";
//			for(int i=0;i<mNumParameters;++i)
//				std::cout<<r[i]<<" ";
//			std::cout<<"\n";

#ifdef ENABLE_EIGEN
			pr = p*r;
			//upate for start solution
			d = pr;
#else
//			pr=r;
			for(int i=0;i<mNumParameters;++i)
			{
//				pr[i]=r[i];
 				pr[i]=p[i]*r[i];
				//upate for start solution
 				d[i]=pr[i];
			}
#endif

#ifdef ENABLE_EIGEN
			alphaNumerator = r.dot(pr);
#else
			alphaNumerator=0;
			for(int i=0;i<mNumParameters;++i)
 				alphaNumerator+=r[i]*pr[i];
#endif
			if (mVerboseLevel>2)
				std::cout<<__FILE__<<" "<<__LINE__<<" normGradient "<<alphaNumerator << " accuracy " <<mAccuracyGradientScaled<<std::endl;

			if (mVerboseLevel>5 && curCycle>0)
				std::cout<< "   Restart after " <<curCycle << " cycles " << std::endl;
			curCycle = 0;
		std::cout.precision(precision);
//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" a1 "<<alphaNumerator<<"\n";
		}


		//begin with cycles bigger 0

		if (mVerboseLevel>2)
			std::cout<<__FILE__<<" "<<__LINE__<<" calc search direction"<<std::endl;


		for(int i=0;i<mNumParameters;++i)
			h[i]=0;

		if(!mMatrixFreeMethod)
			CalculateMatrixVectorProductEBE(d,h);
		else
			CalculateMatrixVectorProductNBN(d,h);
#ifdef ENABLE_EIGEN
	    alphaDenominator = d.dot(h);
		alpha = alphaNumerator/alphaDenominator;
		//set new parameter
		u+= alpha*d;
		r-=alpha*h;
		pr = p*r;
		//no preconditioning
//		pr = r;
		betaNumerator = r.dot(pr);
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;
		d *=beta;
		d +=pr;
#else
	    alphaDenominator=0;
	    squaredNorm=0;
		for(int i=0;i<mNumParameters;++i)
			alphaDenominator+=d[i]*h[i];
		alpha = alphaNumerator/alphaDenominator;
		betaNumerator =0;
		for(int i=0;i<mNumParameters;++i)
		{
				u[i]+= alpha*d[i];
				r[i]-=alpha*h[i];
				pr[i]=p[i]*r[i];
				squaredNorm+=r[i]*r[i];
				betaNumerator+=r[i]*pr[i];
		}
		beta = betaNumerator/ alphaNumerator;
		alphaNumerator = betaNumerator;
		for(int i=0;i<mNumParameters;++i)
		{
			d[i] *=beta;
			d[i] +=pr[i];
		}
#endif
		std::cout.precision(precision);
//		std::cout <<std::setw(width)<<"It.: "<<curIteration<<" i=3 -> h "<<h[3]<<" alpha "<<alpha<<" r "<<r[3]<<" p "<<p[3]<<" pr "<<pr[3]<<" b1 "<<betaNumerator<<" d "<<d[3]<<" u "<<u[3]<<"\n";


		if (mVerboseLevel>2)
		{
			std::cout.precision(precision);
//			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm gradient squared "<<squaredNorm<<std::endl;
			std::cout <<std::setw(width)<<" It.: "<< curIteration<< " norm gradient squared "<<betaNumerator<<std::endl;
			std::cout << std::setw(width)<< "pr " ;
			for (int count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< pr[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "r " ;
			for (int count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< r[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "h " ;
			for (int count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< h[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "d " ;
			for (int count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<< d[count] << "   " ;
			}
			std::cout << std::endl;
			std::cout << std::setw(width)<< "u " ;
			for (int count=0; count<mNumParameters; count++)
			{
				std::cout << std::setw(width)<<u[count] << "   " ;
			}
			std::cout << std::endl;

			std::cout << std::setw(width)<< "alpha "<<alpha<< "beta "<<beta << std::endl;
		}

//		if (squaredNorm<mAccuracyGradientScaled*mAccuracyGradientScaled)
		if (betaNumerator<mAccuracyGradientScaled*mAccuracyGradientScaled)
		{
			if (mVerboseLevel>2)
			{
				std::cout<< "CONVERGED " << std::endl;
			}
			converged = true;
			returnValue = NORMGRADIENT;
			break;
		}

		if (beta<0)
		{
			std::cout<< "Set beta ("<< beta <<") to zero " << std::endl;
			beta=0;
		}


		if (mVerboseLevel>1 && curIteration%mShowSteps==0)
			std::cout<< "Iteration " << curIteration <<" with norm grad squared" << betaNumerator << std::endl;

		//increase iteration and curCycle
		curCycle++;
		curIteration++;

		if (curIteration>mMaxIterations)
		{
			converged = true;
			returnValue = MAXITERATIONS;
			break;
		}

		numFunctionCalls++;
		if (numFunctionCalls>mMaxFunctionCalls)
		{
			converged = true;
			returnValue = MAXFUNCTIONCALLS;
			break;
		}
	}
	isBuild = true;
	// Get force vector
//	r.clear();
//	CalculateReactionForcesEBE(u,r);
//	std::cout<<" Reaction forces \n";
//	for (int count=0; count<mNumParameters; count++)
//	{
//		std::cout << std::setw(width)<< r[count] << "\n   " ;
//	}
//	std::cout<<"\n";

#ifdef SHOW_TIME
    endOpt=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::Optimize] " << difftime(endOpt,startOpt)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
	if (mVerboseLevel>0)
	{
		std::cout<< " "  << std::endl;
		std::cout<< "Number of Function Calls......... " << numFunctionCalls << std::endl;
		std::cout<< "Number of Gradient Calls......... " << numGradientCalls << std::endl;
		std::cout<< "Number of Hessian Calls.......... " << numHessianCalls << std::endl;
		std::cout<< "Number of Iterations............. " << curIteration << std::endl;
		std::cout<< "Active convergence criterion..... " ;
		switch (returnValue)
		{
			case MAXFUNCTIONCALLS:
				std::cout<< "Maximum number of function calls reached." << std::endl;
				break;
			case MAXGRADIENTCALLS:
				std::cout<< "Maximum number of gradient calls reached." << std::endl;
				break;
			case MAXHESSIANCALLS:
				std::cout<< "Maximum number of hessian calls reached." << std::endl;
				break;
			case MAXITERATIONS:
				std::cout<< "Maximum number of iterations reached." << std::endl;
				break;
			case NORMGRADIENT:
				std::cout<< "Norm of preconditioned gradient smaller than prescribed value." << std::endl;
				break;
			default:
				std::cout<< "Unknown convergence criterion." << std::endl;
				break;
		}
		std::cout << std::endl;
		std::cout.precision(precision);
		std::cout << std::setw(width)<< "displacements " ;
		for (int count=0; count<mNumParameters; count++)
		{
			std::cout << std::setw(width)<<u[count] << "   " ;

		}
		std::cout << std::endl;
	}
//	mParameters=u;

	return returnValue;
}

void NuTo::ConjugateGradientGrid::CalcScalingFactors(int& numHessianCalls,myType &p)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif
    //diagonal scaling with scaling factor
	++numHessianCalls;
    double scalefactor=0.00001;
    if(mVerboseLevel>0)
    	std::cout<<"[ConjugateGradientGrid] scale factor "<<scalefactor<<"\n";
    for (int count=0; count<mNumParameters; ++count)
        p[count] *=scalefactor;

//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGrid::CalcScalingFactors] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
//#endif
}

void NuTo::ConjugateGradientGrid::HessianDiag(myType& rHessianDiag)const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    if (mVerboseLevel>2)
    	std::cout<<__FILE__<<" "<<__LINE__<<" in Routine ConjugateGradientGrid::HessianDiag"<<std::endl;
	int numElems=mElemExist.size();

	// global external force vector (active dofs)
    //	myType  force(mNumParameters,1);
	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;elementNumber++)
    {
		if (mElemExist[elementNumber])
		{
			//get the number of the material
			int numStiff = mMaterialOfElem[elementNumber];
			for(int node=0;node<8;++node)
			{
				rHessianDiag[3*mNodeIds[8*elementNumber+node]]+=(mYoungsModulus[numStiff]*mBaseStiffness[0]);
				rHessianDiag[3*mNodeIds[8*elementNumber+node]+1]+=(mYoungsModulus[numStiff]*mBaseStiffness[0]);
				rHessianDiag[3*mNodeIds[8*elementNumber+node]+2]+=(mYoungsModulus[numStiff]*mBaseStiffness[0]);
			}
		}
	}
	for (int i=0;i<mNumParameters;++i)
	{
		if (rHessianDiag[i]!=0)
			rHessianDiag[i]=1./rHessianDiag[i];
	}
//
//#ifdef ENABLE_EIGEN
//	std::cout<<"Hessian: "<<rHessianDiag.transpose()<<std::endl;
//#else
//	std::cout<<"Hessian: ";
//	for(int i=0;i<mNumParameters;++i)
//		std::cout<<rHessianDiag[i]<<" ";
//	std::cout<<std::endl;
//#endif
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::HessianDiag] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
#endif
}

//! @brief ... calculate matrix-vector-product in element-by-element way
//! @brief ... with local vectors
//! @param ... u - prarmeters input, r - gradient output
void NuTo::ConjugateGradientGrid::CalculateMatrixVectorProductEBE(myType &u,myType &r)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME

	int numElems=(int) mElemExist.size();
	int dofsElem=24;
	// global external force vector (active dofs)
    //	myType  force(mNumParameters,1);
	int nodeId=0;
	int numStiff=0;

	std::unique_ptr<double []> residual (new double [24]);
	std::unique_ptr<double []> displacement(new double [24]);

	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;++elementNumber)
    {
		if (mElemExist[elementNumber])
		{
			//get the number of the material
			numStiff = mMaterialOfElem[elementNumber];
			//calculate local return vector with all dofs: r=Ku
			for(int node=0;node<8;++node)
			{
				nodeId=mNodeIds[8*elementNumber+node];
				displacement[3*node]=u[3*nodeId];
				displacement[3*node+1]=u[3*nodeId+1];
				displacement[3*node+2]=u[3*nodeId+2];
			}
			for(int i=0;i<dofsElem;++i)
			{
				residual[i]=0;
				for(int j=0;j<dofsElem;++j)
					residual[i]+=mYoungsModulus[numStiff]*mBaseStiffness[dofsElem*i+j]*displacement[j];
			}
			for(int node=0;node<8;++node)
			{
				nodeId=mNodeIds[8*elementNumber+node];
				if (!mDofIsConstraint[3*nodeId])
					r[3*nodeId]+=residual[3*node];
				if (!mDofIsConstraint[3*nodeId+1])
					r[3*nodeId+1]+=residual[3*node+1];
				if (!mDofIsConstraint[3*nodeId+2])
					r[3*nodeId+2]+=residual[3*node+2];
//			std::cout<<"elem "<<elementNumber<<"  r["<<node<<"] = "<<r[3*nodeId] <<" "<<r[3*nodeId+1] <<" "<<r[3*nodeId+2] <<"\n";
			}
			// Variants-Test:
			// if-clauses after calculation
			// local u and r faster, because of cache
		}
//		else
//			std::cout<<"[CalculateMatrixVectorProductEBE] : "<<elementNumber<<"th element does not exist.\n";
	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateMatrixVectorProductEBE] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

void NuTo::ConjugateGradientGrid::CalculateReactionForcesEBE(myType &u,myType &f)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME

	int numElems=(int) mElemExist.size();
	int dofsElem=24;
	// global external force vector (active dofs)
    //	myType  force(mNumParameters,1);
	int nodeId=0;
	int numStiff=0;

	std::unique_ptr<double []> residual (new double [24]);
	std::unique_ptr<double []> displacement(new double [24]);

	//loop over all elements
	for (int elementNumber=0;elementNumber<numElems;++elementNumber)
    {
		if (mElemExist[elementNumber])
		{
			//get the number of the material
			numStiff = mMaterialOfElem[elementNumber];
			//calculate local return vector with all dofs: f=Ku
			for(int node=0;node<8;++node)
			{
				nodeId=mNodeIds[8*elementNumber+node];
				displacement[3*node]=u[3*nodeId];
				displacement[3*node+1]=u[3*nodeId+1];
				displacement[3*node+2]=u[3*nodeId+2];
			}
			for(int i=0;i<dofsElem;++i)
			{
				residual[i]=0;
				for(int j=0;j<dofsElem;++j)
					residual[i]+=mYoungsModulus[numStiff]*mBaseStiffness[dofsElem*i+j]*displacement[j];
			}
			for(int node=0;node<8;++node)
			{
				nodeId=mNodeIds[8*elementNumber+node];
				f[3*nodeId]+=residual[3*node];
				f[3*nodeId+1]+=residual[3*node+1];
				f[3*nodeId+2]+=residual[3*node+2];
//			std::cout<<"elem "<<elementNumber<<"  f["<<node<<"] = "<<f[3*nodeId] <<" "<<f[3*nodeId+1] <<" "<<f[3*nodeId+2] <<"\n";
			}
			// Variants-Test:
			// if-clauses after calculation
			// local u and f faster, because of cache
		}
	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateMatrixVectorProductEBE] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief ... calculate matrix-vector-product in node-by-node way
//! @brief ... with local vectors
//! @param ... u - prarmeters input, r - gradient output
void NuTo::ConjugateGradientGrid::CalculateMatrixVectorProductNBN(myType &u,myType &r)
{
//#ifdef SHOW_TIME
//    std::clock_t start,end;
//    start=clock();
//#endif //SHOW_TIME
//	std::cout<<"[ConjugatGradientGrid::CalculateMatrixVectorProductNBN]\n";
	int numNodes=(int) mNodeExist.size();
	//loop over all elements
	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
    {
		if (mNodeExist[nodeNumber])
		{
			r[3*nodeNumber+0]=0;
			r[3*nodeNumber+1]=0;
			r[3*nodeNumber+2]=0;
			//27 neighbor nodes
			//calculate local return vector: r=Ku
			for(int i=0;i<27;++i)
			{
				if (mNeighborNodes[27*nodeNumber+i]>-1)
				{
					if (!mDofIsConstraint[3*nodeNumber])
					{
						r[3*nodeNumber+0]+=mEdgeStiffness[27*9*nodeNumber+9*i]*u[3*mNeighborNodes[27*nodeNumber+i]];
						r[3*nodeNumber+0]+=mEdgeStiffness[27*9*nodeNumber+9*i+1]*u[3*mNeighborNodes[27*nodeNumber+i]+1];
						r[3*nodeNumber+0]+=mEdgeStiffness[27*9*nodeNumber+9*i+2]*u[3*mNeighborNodes[27*nodeNumber+i]+2];
					}
					if (!mDofIsConstraint[3*nodeNumber+1])
					{
						r[3*nodeNumber+1]+=mEdgeStiffness[27*9*nodeNumber+9*i+3]*u[3*mNeighborNodes[27*nodeNumber+i]];
						r[3*nodeNumber+1]+=mEdgeStiffness[27*9*nodeNumber+9*i+4]*u[3*mNeighborNodes[27*nodeNumber+i]+1];
						r[3*nodeNumber+1]+=mEdgeStiffness[27*9*nodeNumber+9*i+5]*u[3*mNeighborNodes[27*nodeNumber+i]+2];
					}
					if (!mDofIsConstraint[3*nodeNumber+2])
					{
						r[3*nodeNumber+2]+=mEdgeStiffness[27*9*nodeNumber+9*i+6]*u[3*mNeighborNodes[27*nodeNumber+i]];
						r[3*nodeNumber+2]+=mEdgeStiffness[27*9*nodeNumber+9*i+7]*u[3*mNeighborNodes[27*nodeNumber+i]+1];
						r[3*nodeNumber+2]+=mEdgeStiffness[27*9*nodeNumber+9*i+8]*u[3*mNeighborNodes[27*nodeNumber+i]+2];
					}
//					std::cout<<"node "<<nodeNumber<<" i "<<i<<" neigh "<< mNeighborNodes[27*nodeNumber+i]<<" index stiff "<<27*9*nodeNumber+9*i<<" r["<<nodeNumber<<"] = "<<r[3*nodeNumber] <<" "<<r[3*nodeNumber+1] <<" "<<r[3*nodeNumber+2] <<"\n";
				}
			}
		}
	}
//#ifdef SHOW_TIME
//    end=clock();
//    if (mShowTime)
//        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateMatrixVectorProductNBN] " << difftime(end,start)/CLOCKS_PER_SEC << "sec \n ";
//#endif
}

//! @brief ... calculate start gradient in node-by-node way
//! @brief ... variante II: with 3x3 matrix at node
void NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNode(myType &r)
{
#ifdef ENABLE_MECHANICS
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
	if (mVerboseLevel>2)
		std::cout<<__FILE__<<" "<<__LINE__<<" in Routine CalculateStartGradientNodeByNodeII"<<std::endl;
//	myType gradientNew(mNumParameters,1);
//	int numNodes=mpGrid->GetNumNodes();
//	//nodeIds here ids of all neighbor nodes
//    int* nodeIds;
//	//array with global dof number of each dof of this node
//    int numDofs=3; //for dofs array needed
//
//	// part of matrix one dimensional 9 fields
//	double *array=new double [9*27];
//
//	NodeGrid3D* thisNode=0;
//	// global external force vector (active dofs)
////	myType  force(mpGrid->GetNumDofs(),1);
//	myType locReturn(numDofs,1);
//
//	//loop over all nodes
//	for (int nodeNumber=0;nodeNumber<numNodes;++nodeNumber)
//    {
//		// return vector with all dofs of one node
//		//get pointer to each node
//		thisNode =mpGrid->NodeGridGetNodePtr(nodeNumber);
//		if(thisNode)
//		{
//			locReturn.Resize(numDofs,1);
//			//get belonging nodes
//			nodeIds=thisNode->GetNodeIds();
//	//				std::cout<<" node "<<nodeNumber <<""<<std::endl;
//			for (int node=0;node<27;++node)
//			{
//				//myType* matrix;
//				//get pointer to array of part coefficient matrix for all neighbor nodes
//				array= thisNode->GetPartCoefficient0(0);
//	//			matrix= thisNode->GetPartCoefficientMatrix0(node);
//				//node exists
//		/*
//				if (nodeNumber==62)
//				{
//					std::cout<<" node "<<nodeNumber <<"ids :"<<nodeIds[0]<<"  "<<nodeIds[1]<< std::endl;
//				}
//			*/
//				if (nodeIds[node]>=0)
//				{
//					//myType locDispValues(3,1);
//					//mpGrid->NodeGetDisplacements(nodeIds[node],locDispValues);
//					//calculate nodal return vector and sum up
//					//when global dof is active
//					//when dof is not constraint, then ...
//					//subtract (r=f-Ku) locReturn for active dofs
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+0))
//						gradientNew(nodeNumber*numDofs+0,0)+=array[node*9+0]*r(nodeIds[node]*numDofs,0) + array[node*9+1]*r(nodeIds[node]*numDofs+1,0) + array[node*9+2]*r(nodeIds[node]*numDofs+2,0);
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+1))
//						gradientNew(nodeNumber*numDofs+1,0)+=array[node*9+3]*r(nodeIds[node]*numDofs,0) + array[node*9+4]*r(nodeIds[node]*numDofs+1,0) + array[node*9+5]*r(nodeIds[node]*numDofs+2,0);
//					if (mpGrid->NodeGetConstraintSwitch(nodeNumber*numDofs+2))
//						gradientNew(nodeNumber*numDofs+2,0)+=array[node*9+6]*r(nodeIds[node]*numDofs,0) + array[node*9+7]*r(nodeIds[node]*numDofs+1,0) + array[node*9+8]*r(nodeIds[node]*numDofs+2,0);
//		//				locReturn.operator += (matrix->operator *(locDispValues));
//				/*
//					if (nodeNumber==62)
//					{
//						std::cout<<" node "<<nodeNumber <<" dof 0 locReturn :"<<locReturn(0,0)<<" at loc node "<<nodeIds[node]<< std::endl;
//						std::cout<<" gradOrig :"<<r(nodeIds[node]*numDofs,0)<<" , "<< r(nodeIds[node]*numDofs+1,0)<<" , "<<r(nodeIds[node]*numDofs+2,0)<<std::endl;
//					}
//					*/
//				}
//
//			}
//		}
//    }
//	r=gradientNew;
//		//get global external force vector
//		//! @TODO add load vector
//		// ubdate this routine
//		//mpGrid->BuildGlobalExternalLoadVector(force);
//		//add global external force vector
//		//r+=force;

#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::ConjugateGradientGrid::CalculateStartGradientNodeByNodeII] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
#else
	throw OptimizeException ( "[ConjugateGradientGrid::CalculateStartGradientNodeByNodeII] Modul Optimize is not loaded." );
#endif // ENABLE_MECHANICS
}



#ifdef ENABLE_SERIALIZATION
//! @brief ... save the object to a file
//! @param filename ... filename
//! @param rType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Save ( const std::string &filename, std::string rType)const
{
	try
	{
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
		std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
		std::string tmpStr ( GetTypeId() );
		std::string baseClassStr = tmpStr.substr ( 4,100 );
		if (rType=="BINARY")
		{
			boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
			oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="XML")
		{
			boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
			oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else if (rType=="TEXT")
		{
			boost::archive::text_oarchive ota ( ofs, std::ios::binary );
			ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                & BOOST_SERIALIZATION_NVP(mMaxIterations)
                & BOOST_SERIALIZATION_NVP(mShowSteps);
		}
		else
		{
			throw MathException ( "[ConjugateGradientGrid::Save]File type not implemented." );
		}
	}
	catch ( boost::archive::archive_exception &e )
	{
		std::string s ( std::string ( "[ConjugateGradientGrid::Save]File save exception in boost - " ) +std::string ( e.what() ) );
		std::cout << s << "\n";
		throw MathException ( s );
	}
	catch ( MathException &e )
	{
		throw e;
	}
	catch ( std::exception &e )
	{
		throw MathException ( e.what() );
	}
	catch ( ... )
	{
		throw MathException ( "[Matrix::Save]Unhandled exception." );
	}
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::ConjugateGradientGrid::Restore ( const std::string &filename,  std::string rType)
{
    try
    {
		//transform to uppercase
		std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        std::string tmpString;
		if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oba & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             oxa & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MathException ( "[Matrix::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

            if ( tmpString!=GetTypeId() )
                throw OptimizeException ( "[NuTo::ConjugateGradientGrid::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );

             ota & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Optimizer)
                 & BOOST_SERIALIZATION_NVP(mAccuracyGradient)
                 & BOOST_SERIALIZATION_NVP(mMinDeltaObjBetweenRestarts)
                 & BOOST_SERIALIZATION_NVP(mMaxGradientCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxHessianCalls)
                 & BOOST_SERIALIZATION_NVP(mMaxIterations)
                 & BOOST_SERIALIZATION_NVP(mShowSteps);
        }
		else
		{
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MathException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MathException ( e.what() );
    }
    catch ( ... )
    {
        throw MathException ( "[Matrix::Restore]Unhandled exception." );
    }
}
#endif // ENABLE_SERIALIZATION

//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name ConjugateGradientGrid
std::string NuTo::ConjugateGradientGrid::GetTypeId()const
{
    return std::string("ConjugateGradientGrid");
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ConjugateGradientGrid::Info () const
{
    NuTo::Optimizer::InfoBase();
	std::cout<< "AccuracyGradient" << mAccuracyGradient << std::endl;
	std::cout<< "MinDeltaObjBetweenRestarts" << mMinDeltaObjBetweenRestarts << std::endl;
	std::cout<< "MaxGradientCalls" << mMaxGradientCalls << std::endl;
	std::cout<< "MaxHessianCalls" << mMaxHessianCalls << std::endl;
	std::cout<< "MaxIterations" << mMaxIterations << std::endl;
	std::cout<< "ShowSteps" << mShowSteps << std::endl;
}
void  NuTo::ConjugateGradientGrid::AnsysInput(int rNumParameters,boost::dynamic_bitset<>& elemExist, boost::dynamic_bitset<>& nodeExist,boost::dynamic_bitset<> &rDofIsConstraint,std::vector<double>& youngsModulus,int* rGridDimension,double* rVoxelSpacing,std::vector<int>& materialOfElem,std::vector<int>& allNodesAtVoxel,std::vector<double>& parameters)
{
	// open file
	std::ofstream file;
    file.open("ansysInput");
    file<<"!Ansys Input File \n/prep7 \net,1,solid185 \nkeyopt,1,2,3 \n";

   int countMat=youngsModulus.size();
    for(int i=0;i<countMat;++i)
    {
        file<<"mp,ex,"<<i+1<<","<<youngsModulus[i]<<" \n";
        file<<"mp,prxy,"<<i+1<<",.2 \n";

    }

    for (int z=0;z<rGridDimension[2]+1;++z)
    {
    	for (int y=0;y<rGridDimension[1]+1;++y)
    	{
    		for (int x=0;x<rGridDimension[0]+1;++x)
			{
    			file<<"n,,"<<x*rVoxelSpacing[0]<<","<<y*rVoxelSpacing[1]<<","<<z*rVoxelSpacing[2]<<"\n";
			}
    	}
    }
    int mat=1;
    file<<"type,1 \n mat,"<<mat<<" \n";
    int numElems=(int) elemExist.size();
    int numExistElem(0);
    // attention for more materials
    for (int i=0;i<numElems;++i)
    {
    	if(elemExist[i])
    	{
    		if(materialOfElem[i]!=mat-1)
    		{
    			mat=materialOfElem[i]+1;
    			file<<"mat,"<<mat<<"\n";
    		}
    		file<<"en,"<<i+1<<","<<allNodesAtVoxel[8*i]+1<<","<<allNodesAtVoxel[8*i+1]+1<<","<<allNodesAtVoxel[8*i+2]+1<<","<<allNodesAtVoxel[8*i+3]+1<<","<<allNodesAtVoxel[8*i+4]+1<<","<<allNodesAtVoxel[8*i+5]+1<<","<<allNodesAtVoxel[8*i+6]+1<<","<<allNodesAtVoxel[8*i+7]+1<<"\n";
    		++numExistElem;
    	}
//    	else
//			std::cout<<" AnsysInput: "<<i<<". Element does not exist \n";

    }
    std::cout<<"[NuTo::ConjugateGradientGrid] AnsysInput: "<<numExistElem<<" elements created.\n";
	file<<"nsle,s \n nsel,inve \n ndele,all \n alls \n";
    file<<"alls \n *Get,minx,node,0,mnloc,x \n *Get,miny,node,0,mnloc,y,\n*Get,minz,node,0,mnloc,z,\n";
    file<<" *Get,maxz,node,0,mxloc,z,\n *Get,maxx,node,0,mxloc,x, \n*Get,maxy,node,0,mxloc,y,\n";
	//----------------------------------------------------------------------------------------//
    //Boundary condition: x=0 are constraint
    //Boundary condition: node 0 auy=uz=0 and x=dimx uz=0
	//----------------------------------------------------------------------------------------//
//    file<<"nsel,s,loc,x,minx \n";
//    file<<"d,all,ux,0 \n";
//    file<<"nsel,s,loc,y,miny \n";
//    file<<"nsel,r,loc,z,minz\n";
//    file<<"d,all,uy,0 \n";
//    file<<"d,all,uz,0 \n";
//
//    file<<"nsel,s,loc,x,maxx \n";
//    file<<"nsel,r,loc,y,maxy \n";
//    file<<"nsel,r,loc,z,minz\n";
//    file<<"d,all,uz,0 \n";
//
	//----------------------------------------------------------------------------------------//

    //----------------------------------------------------------------------------------------//
    //Boundary condition: x=max, ux=-1
//	file<<"nsel,s,loc,x,maxx\n";
//	file<<"d,all,ux,-1 \n";
    //----------------------------------------------------------------------------------------//


	//----------------------------------------------------------------------------------------//
    //Boundary condition: all directions for z=0 are constraint
    //Boundary condition:  z=max, uz=-1
	//----------------------------------------------------------------------------------------//    file<<"nsel,s,loc,z,minz\n";
    file<<"nsel,s,loc,z,0\n";
    file<<"d,all,all,0 \n";

    file<<"nsel,s,loc,z,maxz\n";
    // forces
//    file<<"f,all,fz,-1 \n";
    file<<"d,all,uz,-1 \n";
//    file<<"d,all,uz,"<<-(rGridDimension[2]*rVoxelSpacing[2])/20.<<" \n";
	//----------------------------------------------------------------------------------------//    file<<"nsel,s,loc,z,minz\n";
    file.close();
}

