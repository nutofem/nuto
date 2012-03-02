// $Id: Newmark.cpp 575 2011-09-20 18:05:35Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/ptr_container/serialize_ptr_map.hpp>
#endif // ENABLE_SERIALIZATION

# ifdef _OPENMP
#include <omp.h>
# endif

#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseDirectSolverMKLPardiso.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/Newmark.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::Newmark::Newmark ( )  : TimeIntegrationBase ()
{
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::Newmark::Info()const
{
	TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::Newmark::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::Newmark::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::Newmark::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::Newmark::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::Newmark::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::Newmark::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::Newmark::serialize(Archive & ar, const unsigned int version)
{
	#ifdef DEBUG_SERIALIZATION
	    mLogger << "start serialization of Newmark" << "\n";
	#endif
	    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase);

    #ifdef DEBUG_SERIALIZATION
        mLogger << "finish serialization of Newmark" << "\n";
    #endif
}


//! @brief ... save the object to a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Newmark::Save (const std::string &filename, std::string rType )const
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        if(! ofs.is_open())
        {
            throw MechanicsException("[NuTo::Newmark::Save] Error opening file.");
        }

        // write data to file
        std::string typeIdString(this->GetTypeId());
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ("Object_type", typeIdString );
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            std::string tmpString(this->GetTypeId());
            oxa & boost::serialization::make_nvp ("Object_type", typeIdString );
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp("Object_type", typeIdString );
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::Newmark::Save] File type not implemented." );
        }

        // close file
        ofs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::Newmark::Save]File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Newmark::Save] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}


//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::Newmark::Restore (const std::string &filename, std::string rType )
{
#ifdef SHOW_TIME
    std::clock_t start,end;
    start=clock();
#endif
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), (int(*)(int))toupper);

        // open file
        std::ifstream ifs ( filename.c_str(), std::ios_base::binary );
        if(! ifs.is_open())
        {
            throw MechanicsException("[NuTo::Newmark::Restore] Error opening file.");
        }

        std::string typeIdString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Newmark::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oba & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Newmark::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            oxa & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", typeIdString );
            if ( typeIdString != this->GetTypeId() )
            {
                throw MechanicsException ( "[NuTo::Newmark::Restore] Data type of object in file ("+typeIdString+") is not identical to data type of object to read ("+this->GetTypeId() +")." );
            }
            ota & boost::serialization::make_nvp(typeIdString.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[NuTo::Newmark::Restore]File type not implemented" );
        }
        // close file
        ifs.close();
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[NuTo::Newmark::Restore] File save exception in boost - " ) + std::string ( e.what() ) );
        throw MechanicsException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
#ifdef SHOW_TIME
    end=clock();
    if (mShowTime)
        std::cout<<"[NuTo::Newmark::Restore] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << std::endl;
#endif
}

#endif // ENABLE_SERIALIZATION


/*
//! @brief deconstructor
void NuTo::Newmark::Solve(NuTo::StructureBase* rStructure)const
{
    double mBeta(0.25);
    double mGamma(0.5);
    double mConvergenceTolerance(1e-6);
    double mMuM(0.); //Rayleigh damping mass
    double mMuK(0.); //Rayleigh damping stiffness
	double maximumTimeStep(0.01);
	double mMinLineSearchFactor(0.01);
	double mToleranceResidualForce(1e-5);
	double timeStep(0.01);
	double constraintRHS0;

    if (mConstraintLoad!=-1)
    {
    	constraintRHS0 = rStructure->ConstraintGetRHS(mConstraintLoad);
    }
	rStructure->NodeBuildGlobalDofs();

    NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
    rStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
    NuTo::FullMatrix<double> velocitiesActiveDOFsLastConverged(displacementsActiveDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> velocitiesDependentDOFsLastConverged(displacementsDependentDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> accelerationsActiveDOFsLastConverged(displacementsActiveDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> accelerationsDependentDOFsLastConverged(displacementsDependentDOFsLastConverged.GetNumRows(),1);

    //init trial state
    NuTo::FullMatrix<double> displacementsActiveDOFs(displacementsActiveDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> displacementsDependentDOFs(displacementsDependentDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> velocitiesActiveDOFs(velocitiesActiveDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> velocitiesDependentDOFs(velocitiesDependentDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> accelerationsActiveDOFs(accelerationsActiveDOFsLastConverged.GetNumRows(),1);
    NuTo::FullMatrix<double> accelerationsDependentDOFs(accelerationsActiveDOFsLastConverged.GetNumRows(),1);

    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullMatrix<double> dispForceVector;
    NuTo::SparseMatrixCSRVector2General<double> massMatrixCSRVector2;

    //for debugging, you don't need that many vectors
    NuTo::FullMatrix<double> resForceVector, intForceVector, dampForceVector, massForceVector,
                             extForceVector0, extForceVectorDelta, extForceVector,
                             deltaDisplacementsActiveDOFs;

    //allocate solver
    NuTo::SparseDirectSolverMUMPS mySolver;
    //NuTo::SparseDirectSolverMKLPardiso mySolver;

    // build global external load vector at the end of the time step
    rStructure->BuildGlobalExternalLoadVector(extForceVectorDelta);

    //calculate mass matrix (assuming it is not changing during the calculation)
    rStructure->BuildGlobalCoefficientMatrix2(massMatrixCSRVector2);
    std::cout << "mass matrix \n" << NuTo::FullMatrix<double>(massMatrixCSRVector2) << "\n";

    //calculate residual
	rStructure->BuildGlobalGradientInternalPotentialVector(extForceVector0);

	//inertia forces
	massForceVector = massMatrixCSRVector2*accelerationsActiveDOFs;
	extForceVector0 += massForceVector;

	if (mMuM>0)
    {
    	//damping with mass
    	dampForceVector = massMatrixCSRVector2*(velocitiesActiveDOFs * mMuM);
    	extForceVector0 += dampForceVector;
    }
	else
	{
		dampForceVector.Resize(rStructure->GetNumActiveDofs(),1);
	}

    if (mMuK>0)
    {
    	//damping with stiffness
    	rStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
        std::cout << "mass matrix \n" << NuTo::FullMatrix<double>(massMatrixCSRVector2) << "\n";
    	dampForceVector += stiffnessMatrixCSRVector2*(velocitiesActiveDOFs * mMuK);
    	extForceVector0 += dampForceVector;
    }

	extForceVectorDelta -= extForceVector0;

	std::cout << "initialiue the vectors for active acceleration, displacement, and velocities." << "\n";

    double timeCur(mTime0);

    while(timeCur<mTime0+mTimeDelta)
    {
	    //next time step
    	timeCur += timeStep;
	    if (timeCur>mTime0+mTimeDelta)
	    {
	    	timeStep -= mTime0+mTimeDelta;
	    	timeCur = mTime0+mTimeDelta;
	    }
	    if (mConstraintLoad!=-1)
	    {
	    	rStructure->ConstraintSetRHS(mConstraintLoad,constraintRHS0+mConstraintRHSDelta*((timeCur-mTime0)/mTimeDelta));
			//evtl. extract here the dependent dofs of the initial numbering, in order to be able to check for a change in the ordering
		    //since the mass is only build once for all time steps
	    	rStructure->NodeBuildGlobalDofs();
	    }

        //initialize state with assuming a=0 for unconstrained/master dofs (initialize that similar to the static version (dispForce)
        accelerationsActiveDOFs.Resize(accelerationsActiveDOFsLastConverged.GetNumRows(),1);
		velocitiesActiveDOFs = velocitiesActiveDOFsLastConverged + accelerationsActiveDOFsLastConverged * (1.-mGamma);
		displacementsActiveDOFs = displacementsActiveDOFsLastConverged + velocitiesActiveDOFsLastConverged*timeStep+
								  accelerationsActiveDOFsLastConverged * (0.5-mBeta)*timeStep*timeStep;

		//build stiffness (damping is supposed to be a function of the mass matrix and the stiffness)
        rStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);

		//initialize state for constrained dofs
		rStructure->NodeMergeActiveDofValues(displacementsActiveDOFs);
		rStructure->ElementTotalUpdateTmpStaticData();

        //update velocities and accelerations for constrained dofs (maybe not necessary)

		//calculate external force
	    extForceVector = extForceVector0 + extForceVectorDelta*((timeCur-mTime0)/mTimeDelta);

		//internal forces
		rStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
		std::cout << "intForceVector \n" << intForceVector << "\n";
		resForceVector = extForceVector - intForceVector;

		//inertia foces
		massForceVector = massMatrixCSRVector2*accelerationsActiveDOFs;
		resForceVector -= massForceVector;

    	//damping forces
        if (mMuM>0)
        {
        	//damping with mass
        	dampForceVector = massMatrixCSRVector2*(velocitiesActiveDOFs * mMuM);
		    resForceVector -= dampForceVector;
        }
        else
        {
    		dampForceVector.Resize(rStructure->GetNumActiveDofs(),1);
        }

        if (mMuK>0)
        {
        	//damping with stiffness
        	dampForceVector += stiffnessMatrixCSRVector2*(velocitiesActiveDOFs * mMuK);
        	resForceVector -= dampForceVector;
        }

        //norm of residual
        double normPrevResForceVector = resForceVector.Norm();
        double normResForceVector = normPrevResForceVector;

        while (normResForceVector>mConvergenceTolerance)
        {
            //solve the system
            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrixCSRVector2);
            stiffnessMatrixCSR.SetOneBasedIndexing();
            mySolver.Solve(stiffnessMatrixCSR, resForceVector, deltaDisplacementsActiveDOFs);

            //perform a linesearch
            double alpha = 1.;
            double maxResForceVector;
            do
            {
                //add new displacement state
            	NuTo::FullMatrix<double> newDisplacementsActiveDOFs(displacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha);
            	//these vectors are essentially not needed and can be replaced by the direct calculation at the point when they are used
            	NuTo::FullMatrix<double> newVelocitiesActiveDOFs(velocitiesActiveDOFs + deltaDisplacementsActiveDOFs*(alpha*mGamma/(mBeta*timeStep)));
           	    NuTo::FullMatrix<double> newAccelerationsActiveDOFs(accelerationsActiveDOFs + deltaDisplacementsActiveDOFs*(alpha/(mBeta*timeStep*timeStep)));

                //mLogger << " displacementsActiveDOFs" << "\n";
                //displacementsActiveDOFs.Trans().Info(10,3);
                rStructure->NodeMergeActiveDofValues(newDisplacementsActiveDOFs);
                rStructure->ElementTotalUpdateTmpStaticData();

                // calculate residual
        		//calculate residual
        		//internal forces
        		rStructure->BuildGlobalGradientInternalPotentialVector(intForceVector);
        		resForceVector = extForceVector - intForceVector;

        		//inertia foces
        		massForceVector = massMatrixCSRVector2*newAccelerationsActiveDOFs;
        		resForceVector -= massForceVector;

                if (mMuM>0)
                {
                	//damping with mass
                	dampForceVector = massMatrixCSRVector2*(newVelocitiesActiveDOFs * mMuM);
        		    resForceVector -= dampForceVector;
                }

                if (mMuK>0)
                {
                	//damping with stiffness
                	dampForceVector += stiffnessMatrixCSRVector2*(velocitiesActiveDOFs * mMuK);
                	resForceVector -= dampForceVector;
                }

                normResForceVector = resForceVector.Norm();
                maxResForceVector = resForceVector.Abs().Max();

				//PostProcessDataInLineSearch(loadStep, numNewtonIterations, alpha, curLoadFactor, normResidual, normRHS);

                alpha*=0.5;
            }
            while(alpha>mMinLineSearchFactor && normResForceVector>normPrevResForceVector*(1-0.5*alpha) && normResForceVector>mToleranceResidualForce && maxResForceVector>mToleranceResidualForce);

        }

        // update static data
        rStructure->ElementTotalUpdateStaticData();
    }
}
*/

//! @brief solve routine for implicit situations
// todo : check forces due to velocities accelerations on dependent dofs (rhs?) (e.g. for diagonal displacement boundary conditions)
// todo : check for external forces on dependent dofs
NuTo::Error::eError NuTo::Newmark::Solve(NuTo::StructureBase* rStructure)const
{
/*#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    try
    {
        // start analysis
        double timeStep(mMaxDeltaTimeStep);
        double time0 = mTime;

        if (mConstraintLoad!=-1)
        {
        	constraintRHS0 = rStructure->ConstraintGetRHS(mConstraintLoad);
        }
        rStructure->NodeBuildGlobalDofs();

//calculate the initial external load vector from the initial residual, to be done
NuTo::FullMatrix<double> extForceVector0(rStructure->GetNumActiveDofs(),1);
exit(0);

        NuTo::FullMatrix<double> displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged;
        rStructure->NodeExtractDofValues(displacementsActiveDOFsLastConverged,displacementsDependentDOFsLastConverged);
        NuTo::FullMatrix<double> velocitiesActiveDOFsLastConverged,velocitiesDependentDOFsLastConverged;
        rStructure->NodeExtractDofFirstTimeDerivativeValues(velocitiesActiveDOFsLastConverged,velocitiesDependentDOFsLastConverged);
        NuTo::FullMatrix<double> accelerationsActiveDOFsLastConverged,accelerationsDependentDOFsLastConverged;
        rStructure->NodeExtractDofSecondTimeDerivativeValues(accelerationsActiveDOFsLastConverged,accelerationsDependentDOFsLastConverged);

        InitBeforeNewLoadStep(loadStep);
        if (mNumActiveDofs==0)
        {
        	exit(0);
        }

        //init some auxiliary variables
        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
        NuTo::SparseMatrixCSRVector2Symmetric<double> massMatrixCSRVector2;
        NuTo::SparseMatrixCSRVector2Symmetric<double> hessianMatrixCSRVector2;

        NuTo::FullMatrix<double> displacementsActiveDOFs;
        NuTo::FullMatrix<double> accelerationsActiveDOFs;
        NuTo::FullMatrix<double> velocitiesActiveDOFs;
        NuTo::FullMatrix<double> dispIntForceVector;
        NuTo::FullMatrix<double> dispMassForceVector;
        NuTo::FullMatrix<double> intForceVector;
        NuTo::FullMatrix<double> massForceVector;
        NuTo::FullMatrix<double> extForceVector;
        NuTo::FullMatrix<double> rhsVector;

        // build global external load vector and RHS vector this is assumed to be independent of the displacement state
        this->BuildGlobalExternalLoadVector(deltaExtForceVector);
        deltaExtForceVector-=extForceVector0;

        //allocate solver
        NuTo::SparseDirectSolverMUMPS mySolver;
        //NuTo::SparseDirectSolverMKLPardiso mySolver;
#ifdef SHOW_TIME
        if (mShowTime==true)
            mySolver.SetShowTime(true);
        else
            mySolver.SetShowTime(false);
#endif
        bool convergenceInitialLoadStep(false);
        //just try to calculate the stiffness matrix (at t=t0 and the internal potential for the first step t=t0+deltaT)
        while (convergenceInitialLoadStep==false)
        {
            try
            {
                //calculate stiffness
                mTime=time0 + timeStep;
                if (mConstraintLoad!=-1)
                {
        	    	rStructure->ConstraintSetRHS(mConstraintLoad,constraintRHS0+mConstraintRHSDelta*((mTime-time0)/mTimeDelta));
                }
                rStructure->NodeBuildGlobalDofs();

                error = rStructure->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispIntForceVector);
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                    {
                        //decrease load step
                        timeStep*=mDecreaseFactor;

                        //check for minimum delta (this mostly indicates an error in the software
                        if (timeStep<mMinDeltaTimeStep)
                        {
                            mLogger << "[NuTo::StructureBase::NewtonRaphson] No convergence for BuildGlobalCoefficientMatrix0 calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                            return Error::NO_CONVERGENCE;
                        }
                        continue;
                    }
                    else
                        return error;
                }
                error = rStructure->BuildGlobalCoefficientMatrix2(massMatrixCSRVector2, dispMassForceVector);
                if (error!=Error::SUCCESSFUL)
                {
                	throw MechanicsException("[NuTo::TimeIntegration::Solve] error building mass matrix.");
                }

                //since that factor includes the relation between the delta disp and delta acceleration
                dispDampForceVector = dispMassForceVector * (mDampingMass * mGamma/(mBeta*timeStep));

                //since that factor includes the relation between the delta disp and delta acceleration
                dispMassForceVector *= 1./(mBeta*timeStep*timeStep);

                //update displacements of all nodes according to the new conre mat
                {
                    NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                    NuTo::FullMatrix<double> displacementsDependentDOFsCheck1;
                    this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                    this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                    //initialize velocities and accelerations of the current time step based on the previous one
                    accelerationsActiveDOFs = (accelerationsActiveDOFsLastConverged*(0.5-mBeta)*timeStep+velocitiesActiveDOFsLastConverged)*(-1./(timeStep*mBeta));
                    velocitiesActiveDOFs = velocitiesActiveDOFsLastConverged + accelerationsActiveDOFsLastConverged * (1.-mGamma)*timeStep + accelerationsActiveDOFs * (mGamma*timeStep);
                    this->NodeMergeActiveDofFirstTimeDerivativeValues(velocitiesActiveDOFs);
                    this->NodeMergeActiveDofSecondTimeDerivativeValues(accelerationsActiveDOFs);
                    error = this->ElementTotalUpdateTmpStaticData();
                    if (error!=Error::SUCCESSFUL)
                    {
                        if (error==Error::NO_CONVERGENCE)
                        {
                            //decrease time step
                            timeStep*=mDecreaseFactor;

                            //restore initial state
                            if (mConstraintLoad!=-1)
                            {
                    	    	rStructure->ConstraintSetRHS(mConstraintLoad,constraintRHS0);
                            }
                            this->NodeBuildGlobalDofs();
                            this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                            this->NodeMergeActiveDofFirstTimeDerivativeValues(velocitiesActiveDOFsLastConverged);
                            this->NodeMergeActiveDofSecondTimeDerivativeValues(accelerationsActiveDOFsLastConverged);

                            //check for minimum delta (this mostly indicates an error in the software
                            if (timeStep<mMinDeltaTimeStep)
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson] return with no convergence" << "\n";
                                return Error::NO_CONVERGENCE;
                            }
                            continue;
                        }
                        else
                            return error;
                    }
                }

                //calculate the residual for the new state (with applied constraints)
                error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                if (error!=Error::SUCCESSFUL)
                {
                    if (error==Error::NO_CONVERGENCE)
                    {
                        //decrease load step
                        timeStep*=mDecreaseFactor;

                        //restore initial state
                        if (mConstraintLoad!=-1)
                        {
                	    	rStructure->ConstraintSetRHS(mConstraintLoad,constraintRHS0);
                        }
                        this->NodeBuildGlobalDofs();
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                        this->NodeMergeActiveDofFirstTimeDerivativeValues(velocitiesActiveDOFsLastConverged);
                        this->NodeMergeActiveDofSecondTimeDerivativeValues(accelerationsActiveDOFsLastConverged);

                        //check for minimum delta (this mostly indicates an error in the software
                        if (timeStep<mMinDeltaTimeStep)
                        {
                            mLogger << "No convergence for initial resforce/stiffness calculation, delta strain factor for initial increment smaller than minimum " << "\n";
                            return Error::NO_CONVERGENCE;
                        }
                        continue;
                    }
                    else
                        return error;
                }

                dampForceVector = (massMatrixCSRVector2 * velocitiesActiveDOFs) * mDampingMass;
                massForceVector = massMatrixCSRVector2 * accelerationsActiveDOFs;

                convergenceInitialLoadStep = true;
            }
            catch(MechanicsException& e)
            {
                e.AddMessage("[NuTo::StructureBase::NewtonRaphson] Error in Newton-Raphson iteration.");
                throw e;
            }
        }
        extForceVector = extForceVector0 + deltaExtForceVector*(timeStep/mTimeDelta);

        rhsVector = extForceVector - intForceVector - dampForceVector - massForceVector;

        //attention this is only different for the first iteration step
        //since the internal force due to the applied constraints is not considered for the first iteration
        //in order to balance it (no localization in the boundary region)
        //for the linesearch this internal force has to be considered in order to obtain for a linesearch
        //factor of zero the normRHS
        double normRHS = rhsVector.Norm();
        rhsVector = extForceVector + dispIntForceVector + dispDampForceVector + dispMassForceVector;

        //calculate absolute tolerance for matrix entries to be not considered as zero
        double maxValue, minValue, ToleranceZeroStiffness;
        if (stiffnessMatrixCSRVector2.GetNumColumns()==0)
        {
            maxValue = 1.;
            minValue = 1.;
        }
        else
        {
            stiffnessMatrixCSRVector2.Max(maxValue);
            stiffnessMatrixCSRVector2.Min(minValue);
        }
        //mLogger << "min and max " << minValue << " , " << maxValue << "\n";

        ToleranceZeroStiffness = (1e-14) * (fabs(maxValue)>fabs(minValue) ?  fabs(maxValue) : fabs(minValue));
        this->SetToleranceStiffnessEntries(ToleranceZeroStiffness);
        //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
        //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
        //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

        //mySolver.ExportVtkDataFile(std::string("/home/unger3/develop/nuto_build/examples/c++/FineScaleConcurrentMultiscale") + std::string("0") + std::string(".vtk"));

        //store the structure only once in order to be able to restore the situation before entering the routine
        rIsSaved = false;

        //repeat until final time is reached
        bool convergenceStatusLoadSteps(false);
        while (!convergenceStatusLoadSteps)
        {
            double normResidual(1);
            double maxResidual(1);
            int numNewtonIterations(0);
            double alpha(1.);
            TimeIntegration::eConvergenceState convergenceStatus(TimeIntegration::CONTINUE_ITERATION);
            //0 - not converged, continue Newton iteration
            //1 - converged
            //2 - stop iteration, decrease load step
            while(convergenceStatus==TimeIntegration::CONTINUE_ITERATION)
            {
                numNewtonIterations++;

                if (numNewtonIterations>mMaxNumNewtonIterations && alpha<0.25)
                {
                    if (mVerboseLevel>5)
                    {
                        mLogger << "numNewtonIterations (" << numNewtonIterations << ") > MAXNUMNEWTONITERATIONS (" << mMaxNumNewtonIterations << ")" << "\n";
                    }
                    convergenceStatus = TimeIntegration::DECREASE_TIME_STEP;
                    break;
                }

                // solve
                NuTo::FullMatrix<double> deltaDisplacementsActiveDOFs;
                NuTo::SparseMatrixCSRGeneral<double> hessianMatrixCSR(stiffnessMatrixCSRVector2+massMatrixCSRVector2*(1./(mBeta*mtimeStep*mtimeStep)+mGamma*mDampingMass/(mBeta*mTimeStep)));
                stiffnessMatrixCSR.SetOneBasedIndexing();
                try
                {
                    mySolver.Solve(stiffnessMatrixCSR, rhsVector, deltaDisplacementsActiveDOFs);
                }
                catch(...)
                {
                    mLogger << "Error solving system of equations using mumps." << "\n";
                    if (mNumActiveDofs<1000)
                    {
                        NuTo::FullMatrix<double> hessianMatrixFull(hessianMatrixCSR);
                        if (mNumActiveDofs<30)
                        {
                            mLogger << "hessian full" << "\n";
                            mLogger.Out(hessianMatrixFull,12,3);
                        }
                        NuTo::FullMatrix<double> eigenValues;
                        hessianMatrixFull.EigenValuesSymmetric(eigenValues);
                        mLogger << "eigenvalues" << "\n";
                        mLogger.Out(eigenValues.Trans(),12,3);
                        NuTo::FullMatrix<double> eigenVectors;
                        hessianMatrixFull.EigenVectorsSymmetric(eigenVectors);
                        mLogger << "eigenvector 1" << "\n";
                        mLogger.Out(eigenVectors.GetColumn(0).Trans(),12,3);
                    }
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] Error solving system of equations using mumps.");
                }

                //mLogger << " rhsVector" << "\n";
                //rhsVector.Trans().Info(10,3);
                //std::cout << " delta_disp" << "\n";
                //deltaDisplacementsActiveDOFs.Trans().Info(10,3);

                //perform a linesearch
                NuTo::FullMatrix<double> oldDisplacementsActiveDOFs;
                NuTo::FullMatrix<double> oldVelocitiesActiveDOFs;
                NuTo::FullMatrix<double> oldAccelerationsActiveDOFs;
                oldDisplacementsActiveDOFs = displacementsActiveDOFs;
                oldVelocitiesActiveDOFs = velocitiesActiveDOFs;
                oldAccelerationsActiveDOFs = accelerationsActiveDOFs;

                alpha = 1.;
                do
                {
                    //add new displacement state
                    displacementsActiveDOFs = oldDisplacementsActiveDOFs + deltaDisplacementsActiveDOFs*alpha;
                    velocitiesActiveDOFs    = oldVelocitiesActiveDOFs    + deltaDisplacementsActiveDOFs*(alpha*mGamma/(mBeta*timeStep));
                    accelerationsActiveDOFs = oldAccelerationsActiveDOFs + deltaDisplacementsActiveDOFs*(alpha/(mBeta*timeStep*timeStep));

                    //mLogger << " displacementsActiveDOFs" << "\n";
                    //displacementsActiveDOFs.Trans().Info(10,3);
                    this->NodeMergeActiveDofValues(displacementsActiveDOFs);
//this merge is evtl. not necessary
                    this->NodeMergeActiveDofFirstTimeDerivativeValues(velocitiesActiveDOFs);
                    this->NodeMergeActiveDofSecondTimeDerivativeValues(accelerationsActiveDOFs);
                    try
                    {
                        Error::eError error = this->ElementTotalUpdateTmpStaticData();
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                convergenceStatus=TimeIntegration::DECREASE_TIME_STEP;
                                mLogger << "Constitutive model is not converging, try with smaller load step" << "\n";
                            }
                            else
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson] error while calling the gradient (resforce) routine." << "\n";
                                return error;
                            }
                        }
                        else
                        {
                        	// calculate residual
                            Error::eError error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                            if (error!=Error::SUCCESSFUL)
                            {
                                if (error==Error::NO_CONVERGENCE)
                                {
                                    convergenceStatus=TimeIntegration::DECREASE_TIME_STEP;
                                    mLogger << "Constitutive model is not converging, try with smaller load step" << "\n";
                                }
                                else
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson] error while calling the gradient (resforce) routine." << "\n";
                                    return error;
                                }
                            }
                            else
                            {
                                //mLogger << "intForceVector "  << "\n";
                                //intForceVector.Trans().Info(10,3);
                                rhsVector = extForceVector - intForceVector;

								//damping with mass
								rhsVector -= massMatrixCSRVector2*(velocitiesActiveDOFs * mDampingMass);

								//inertia forces
                               	rhsVector -= massMatrixCSRVector2*accelerationsActiveDOFs;

                               	normResidual = rhsVector.Norm();
                                maxResidual = rhsVector.Abs().Max();

                                //mLogger << "total energy of system " << ElementTotalGetTotalEnergy() << "\n";
                                PostProcessDataInLineSearch(loadStep, numNewtonIterations, alpha, curLoadFactor, normResidual, normRHS);
                            }
                        }
                    }
                    catch(MechanicsException& e)
                    {
                        e.AddMessage("[NuTo::StructureBase::NewtonRaphson] NuTo::MechanicsException while calling the gradient (resforce) routine.");
                        throw e;
                    }
                    catch(...)
                    {
                        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] Error calling the gradient (resforce) routine.");
                    }

                    alpha*=0.5;

                    if (alpha<mMinLineSearchFactor)
                    	convergenceStatus=NuTo::TimeIntegration::DECREASE_TIME_STEP;
                }
                while(convergenceStatus!=NuTo::TimeIntegration::DECREASE_TIME_STEP && normResidual>normRHS*(1-0.5*alpha) && normResidual>mToleranceResidualForce && maxResidual>mToleranceResidualForce);

                rStructure->PostProcessDataAfterLineSearch(loadStep, numNewtonIterations, 2.*alpha, curLoadFactor, normResidual, rhsVector);

                if (convergenceStatus==NuTo::TimeIntegration::DECREASE_TIME_STEP)
                    break;

                //mLogger << "\n" << "Newton iteration " << numNewtonIterations << ", final alpha " << 2*alpha << ", normResidual " << normResidual<< ", maxResidual " << maxResidual<<"\n";

                //check convergence
                if (normResidual<mToleranceResidualForce || maxResidual<mToleranceResidualForce)
                {
                    //this test is only relevant for problems with a model adaptation, otherwise, just assume a converged solution and continue
                    if(rStructure->AdaptModel())
                    {
                        throw MechanicsException("[NuTo::Newmark::Solve] adaptation not yet implemented.");
                        //change ext force vector to get equilibrium
                    }
                    else
                    {
                        rStructure->PostProcessDataAfterConvergence(loadStep, numNewtonIterations, mTime, timeStep, normResidual);
                        convergenceStatus=TimeIntegration::CONVERGED;
                        //CheckStiffness();
                        //NodeInfo(12);
                        break;
                    }
                }

                //convergence status == 0 (continue Newton iteration)
                normRHS = rhsVector.Norm();
                //build new stiffness matrix
                //there should not be a problem here, since the gradient has already been successfully caclulated for that loading scenario
                Error::eError error = this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
                if (error!=Error::SUCCESSFUL)
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] stiffness for first step of new load increment could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");

                //mLogger << dispForceVector.Norm() << "\n";
//check stiffness
//CheckStiffness();
                //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
                //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
                //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";
            }

            if (convergenceStatus==TimeIntegration::CONVERGED)
            {
                // the update is only required to allow for a stepwise solution procedure in the fine scale model
                // a final update is only required for an update on the macroscale, otherwise,the original state has
                // to be reconstructed.

                if (mTime>time0+mTimeDelta-1-1e-8)
                {
                    if (rSaveStructureBeforeUpdate==false)
                    {
                        Error::eError error = this->ElementTotalUpdateStaticData();
                        if (error!=Error::SUCCESSFUL)
                            throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] update could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");
                        this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                    }
                    convergenceStatusLoadSteps=true;
                }
                else
                {
                    Error::eError error = this->ElementTotalUpdateStaticData();
                    if (error!=Error::SUCCESSFUL)
                        throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] update could not be calculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");

                    this->PostProcessDataAfterUpdate(loadStep, numNewtonIterations, curLoadFactor, deltaLoadFactor, normResidual);
                    //store the last converged step in order to be able to go back to that state
                    displacementsActiveDOFsLastConverged  = displacementsActiveDOFs;

                    //eventually increase load step
                    if (mAutomaticLoadstepControl)
                    {
                        if (numNewtonIterations<mMinNumNewtonIterations)
                        {
                            timeStep*=mIncreaseFactor;
                        }
                        if (timeStep>mMaxTimeStep)
                        	timeStep = mMaxTimeStep;
                    }

                    //increase displacement
                    mTime+=timeStep;
                    if (mTime>time0+mTimeDelta)
                    {
                    	timeStep -= mTime - (time0+mTimeDelta);
                    	mTime= time0+mTimeDelta;
                    }
                }
                loadStep++;
                //initialize some stuff before a new load step (e.g. output directories for visualization, if required)
                InitBeforeNewLoadStep(loadStep);
            }
            else
            {
                assert(convergenceStatus==TimeIntegration::DECREASE_TIME_STEP);
                if (mAutomaticLoadstepControl==false)
                    return Error::NO_CONVERGENCE;

                mLogger << "no convergence with current time step (" << timeStep << "), current not converging load factor " << mTime << "\n";
                //mLogger << "check stiffness " << "\n";
                //CheckStiffness();
                mLogger << "and continue with smaller time step " << "\n";

                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                mTime-=timeStep;

                //restore previous converged state
                if (mConstraintLoad!=-1)
                {
        	    	rStructure->ConstraintSetRHS(mConstraintLoad,constraintRHS0+mConstraintRHSDelta*((mTime-time0)/mTimeDelta));
                }

                // build global dof numbering
                this->NodeBuildGlobalDofs();

                //set previous converged displacements
                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
//this merge is evtl. not necessary
				this->NodeMergeActiveDofFirstTimeDerivativeValues(velocitiesActiveDOFsLastConverged);
				this->NodeMergeActiveDofSecondTimeDerivativeValues(accelerationsActiveDOFsLastConverged);
                 Error::eError error = this->ElementTotalUpdateTmpStaticData();
                if (error!=Error::SUCCESSFUL)
                    throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] last converged state could not be recalculated, but since the gradient at that point has been evaluated successfully, there is a problem in the implementation.");
                //this first part of the routine is only relevant for the multiscale model, since an update on the fine scale should only be performed
                //for an update on the coarse scale
                //as a consequence, in an iterative solution with updates in between, the initial state has to be restored after leaving the routine
                if (rSaveStructureBeforeUpdate==true && rIsSaved==false)
                {
                    assert(mTime==time0);
                    //store the structure only once in order to be able to restore the situation before entering the routine
                    this->SaveStructure(rSaveStringStream);
                    rIsSaved = true;
                }

                //decrease time step
                timeStep*=mDecreaseFactor;
                mTime+=timeStep;

                //check for minimum delta (this mostly indicates an error in the software)
                if (timeStep<mMinDeltaTimeStep)
                {
                    mLogger << "return with a MechanicsNoConvergenceException " << "\n";
                    return Error::NO_CONVERGENCE;
                }
            }

            if (!convergenceStatusLoadSteps)
            {
                try
                {
                    bool convergenceConstitutive(false);
                    while (convergenceConstitutive==false)
                    {
                        this->SetLoadFactor(curLoadFactor);
                        this->NodeBuildGlobalDofs();

                        //update stiffness in order to calculate new dispForceVector, but still with previous displacement state
                        Error::eError error = this->BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (mAutomaticLoadstepControl==false)
                                return Error::NO_CONVERGENCE;
                            mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                            //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                            //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                            //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                            curLoadFactor-=deltaLoadFactor;

                            //set the previous displacement state
                            this->SetLoadFactor(curLoadFactor);

                            // build global dof numbering
                            this->NodeBuildGlobalDofs();

                            //set previous converged displacements
                            this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                            Error::eError error = this->ElementTotalUpdateTmpStaticData();
                            if (error!=Error::SUCCESSFUL)
                                throw MechanicsException("[NuTo::StructureBase::NewtonRaphson] previous converged load step could not be recalculated, check your implementation.");

                            //decrease load step
                            deltaLoadFactor*=mDecreaseFactor;
                            curLoadFactor+=deltaLoadFactor;

                            //check for minimum delta (this mostly indicates an error in the software
                            if (deltaLoadFactor<mMinDeltaLoadFactor)
                            {
                                mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                return Error::NO_CONVERGENCE;
                            }
                            continue;
                        }
                        //int numRemoved = stiffnessMatrixCSRVector2.RemoveZeroEntries(ToleranceZeroStiffness,0);
                        //int numEntries = stiffnessMatrixCSRVector2.GetNumEntries();
                        //mLogger << "stiffnessMatrix: num zero removed " << numRemoved << ", numEntries " << numEntries << "\n";

                        //update displacements of all nodes according to the new conre mat
                        NuTo::FullMatrix<double> displacementsActiveDOFsCheck;
                        NuTo::FullMatrix<double> displacementsDependentDOFsCheck;
                        this->NodeExtractDofValues(displacementsActiveDOFsCheck, displacementsDependentDOFsCheck);
                        this->NodeMergeActiveDofValues(displacementsActiveDOFsCheck);
                        error = this->ElementTotalUpdateTmpStaticData();
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                if (mAutomaticLoadstepControl==false)
                                    return Error::NO_CONVERGENCE;
                                mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                                curLoadFactor-=deltaLoadFactor;

                                //set the previous displacement state
                                this->SetLoadFactor(curLoadFactor);

                                // build global dof numbering
                                this->NodeBuildGlobalDofs();

                                //set previous converged displacements
                                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                                if (error!=Error::SUCCESSFUL)
                                    return error;

                                //decrease load step
                                deltaLoadFactor*=mDecreaseFactor;
                                curLoadFactor+=deltaLoadFactor;

                                //check for minimum delta (this mostly indicates an error in the software
                                if (deltaLoadFactor<mMinDeltaLoadFactor)
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                    mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                    return Error::NO_CONVERGENCE;
                                }
                                continue;
                            }
                            else
                                return error;
                        }

                        // calculate initial residual for next load step
                        //mLogger<<" calculate gradient 1532" << "\n";
                        error = this->BuildGlobalGradientInternalPotentialVector(intForceVector);
                        if (error!=Error::SUCCESSFUL)
                        {
                            if (error==Error::NO_CONVERGENCE)
                            {
                                if (mAutomaticLoadstepControl==false)
                                    return Error::NO_CONVERGENCE;
                                mLogger << "Constitutive model is not converging in initial try after load increment, I'll try with smaller load step" << "\n";

                                //calculate stiffness of previous loadstep (used as initial stiffness in the next load step)
                                //this is done within the loop in order to ensure, that for the first step the stiffness matrix of the previous step is used
                                //otherwise, the additional boundary displacements will result in an artifical localization in elements at the boundary
                                curLoadFactor-=deltaLoadFactor;

                                //set the previous displacement state
                                this->SetLoadFactor(curLoadFactor);

                                // build global dof numbering
                                this->NodeBuildGlobalDofs();

                                //set previous converged displacements
                                this->NodeMergeActiveDofValues(displacementsActiveDOFsLastConverged);
                                Error::eError error = this->ElementTotalUpdateTmpStaticData();
                                if (error!=Error::SUCCESSFUL)
                                    return error;

                                //decrease load step
                                deltaLoadFactor*=mDecreaseFactor;
                                curLoadFactor+=deltaLoadFactor;

                                //check for minimum delta (this mostly indicates an error in the software
                                if (deltaLoadFactor<mMinDeltaLoadFactor)
                                {
                                    mLogger << "[NuTo::StructureBase::NewtonRaphson]: No convergence after the initial stiffness calculation after a load decrease.";
                                    mLogger << "curLoadFactor " << curLoadFactor << ", deltaLoadFactor " << deltaLoadFactor << "mMinDeltaLoadFactor" << mMinDeltaLoadFactor<< "\n";
                                    return Error::NO_CONVERGENCE;
                                }
                                continue;
                            }
                            else
                                return error;
                        }

                        convergenceConstitutive = true;
                    }

                }
                catch(NuTo::MechanicsException& e)
                {
                    e.AddMessage("[NuTo::StructureBase::NewtonRaphson]: error calculating new displacement increment.");
                    throw e;
                }
                catch(...)
                {
                    throw NuTo::MechanicsException("[NuTo::StructureBase::NewtonRaphson]: Exception after the initial stiffness calculation after a load decrease.");
                }

                //update rhs vector for next Newton iteration
                rhsVector = extForceVector - intForceVector;
                normRHS = rhsVector.Norm();
                //attention this is only different for the first iteration step (load application)
                //since the internal force due to the applied constraints is not considered for the first iteration
                //in order to balance it (no localization in the boundary region)
                //for the linesearch this internal force has to be considered in order to obtain for a linesearch
                //factor of zero the normRHS
                rhsVector = dispForceVector + extForceVector;
                if (rInitialStateInEquilibrium==false && curLoadFactor==deltaLoadFactor)
                {
                    rhsVector -= intForceVectorInit;
                }
            }
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::Newmark::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        rStructure->GetLogger()<<"[NuTo::Newmark::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
    	rStructure->GetLogger()<< "[NuTo::Newmark::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
*/
    return NuTo::Error::SUCCESSFUL;

}


#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::Newmark)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
