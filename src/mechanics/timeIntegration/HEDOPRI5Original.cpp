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

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/HEDOPRI5Original.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"



NuTo::HEDOPRI5Original::HEDOPRI5Original (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mTimeStep = 0.;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::HEDOPRI5Original::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::HEDOPRI5Original::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
	return 2.8/std::sqrt(maxGlobalEigenValue);
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
double NuTo::HEDOPRI5Original::GetStageTimeFactor(int rStage)const
{
    assert(rStage<8);
	double s;

	switch(rStage)
	{
	case 0:
		s = 0.;
		break;
	case 1:
        s = 1./5.;
		break;
	case 2:
        s = 3./10.;
		break;
	case 3:
        s = 4./5.;
		break;
    case 4:
        s = 8./9.;
        break;
    case 5:
        s = 1.0;
        break;
    case 6:
        s = 1.0;
        break;
    case 7:
        s = 19./20.;
        break;
	default:
        throw MechanicsException ( "[NuTo::HEDOPRI5Original::GetStageTimeFactor] rStage>7 not implemented." );
	}
	return s;
}

//! @brief ... return delta time factor of intermediate stages (c in Butcher tableau, but only the delta to the previous step)
// so essentially it's c_n-c_(n-1)
bool NuTo::HEDOPRI5Original::HasTimeChanged(int rStage)const
{
    assert(rStage<8);
	bool s;
	switch(rStage)
	{
	case 0:
		s = false; //same as last step from the last iteration
		break;
	case 1:
		s = true;
		break;
	case 2:
        s = true;
		break;
	case 3:
        s = true;
		break;
    case 4:
        s = true;
        break;
    case 5:
        s = true;
        break;
    case 6:
        s = true;
        break;
    case 7:
        s = true;
        break;
	default:
        throw MechanicsException ( "[NuTo::HEDOPRI5Original::HasTimeChanged] rStage>7 not implemented." );
	}
	return s;
}


//! @brief ... return scaling for the intermediate stage for y (a in Butcher tableau)
void NuTo::HEDOPRI5Original::GetStageDerivativeFactor(std::vector<double>& rWeight, int rStage)const
{
    assert(rStage<8);
    assert(rWeight.size()==7);
	switch(rStage)
	{
	case 0:
		break;
	case 1:
        rWeight[0] = 1./5.;
		break;
	case 2:
        rWeight[0] = 3./40.;
        rWeight[1] = 9./40.;
		break;
	case 3:
        rWeight[0] =  44./45.;
        rWeight[1] = -56./15.;
        rWeight[2] =  32./9.;
		break;
    case 4:
        rWeight[0] =  19372./6561.;
        rWeight[1] = -25360./2187.;
        rWeight[2] =  64448./6561.;
        rWeight[3] =  -212./729.;
        break;
    case 5:
        rWeight[0] =  9017./3168.;
        rWeight[1] = -355./33.;
        rWeight[2] =  46732./5247.;
        rWeight[3] =  49./176.;
        rWeight[4] =  -5103./18656.;
        break;
    case 6:
        rWeight[0] =  35./384.;
        rWeight[1] =  0;
        rWeight[2] =  500./1113.;
        rWeight[3] =  125./192.;
        rWeight[4] =  -2187./6784.;
        rWeight[5] =  11./84.;
        break;
    case 7:
        rWeight[0] = -18611506045861./19738176307200.;
        rWeight[1] =  59332529./14479296.;
        rWeight[2] = -2509441598627./893904224850.;
        rWeight[3] =  2763523204159./3289696051200.;
        rWeight[4] = -41262869588913./116235927142400.;
        rWeight[5] =  46310205821./287848404480.;
        rWeight[6] = -3280./75413.;
        break;
	default:
        throw MechanicsException ( "[NuTo::HEDOPRI5Original::GetStageDerivativeFactor] rStage>6 not implemented." );
	}
}

//! @brief ... return weights for the intermediate stage for y (b in Butcher tableau)
double NuTo::HEDOPRI5Original::GetStageWeights(int rStage)const
{
    assert(rStage<8);
	double s;
	switch(rStage)
	{
	case 0:
        s = 35./384.;
		break;
	case 1:
        s = 0.;
		break;
	case 2:
        s = 500./1113.;
		break;
	case 3:
        s = 125./192.;
		break;
    case 4:
        s = -2187./6784.;
        break;
    case 5:
        s = 11./84.;
        break;
    case 6:
        s = 0.;
        break;
    case 7:
        s = 0.;
        break;
	default:
        throw MechanicsException ( "[NuTo::HEDOPRI5Original::GetStageWeights] rStage>7 not implemented." );
	}
	return s;
}
#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::HEDOPRI5Original::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::HEDOPRI5Original::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::HEDOPRI5Original::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::HEDOPRI5Original::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::HEDOPRI5Original::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::HEDOPRI5Original::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::HEDOPRI5Original::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of HEDOPRI5Original" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RungeKuttaBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of HEDOPRI5Original" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the file
//!            in case of restoring from a file with the wrong object type, the file id is printed
//! @return    class name
std::string NuTo::HEDOPRI5Original::GetTypeId()const
{
    return std::string("HEDOPRI5Original");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::RungeKuttaBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::RungeKuttaBase::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of RungeKuttaBase" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mTimeStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of RungeKuttaBase" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::Error::eError NuTo::HEDOPRI5Original::Solve(double rTimeDelta)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime();
#endif
    start=clock();
#endif
    try
    {
        if (mTimeStep==0.)
        {
            if (this->HasCriticalTimeStep())
            {
                mTimeStep = this->CalculateCriticalTimeStep();
            }
            else
            {
                throw MechanicsException("[NuTo::RungeKuttaBase::Solve] time step not set for unconditional stable algorithm.");
            }
        }

        std::cout << "modify computation of critical time step to include the dependence on the time integration scheme." << std::endl;
        //calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time step

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;

        //renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

        //calculate constraint matrix
        NuTo::SparseMatrixCSRGeneral<double> CmatTmp;

        mStructure->ConstraintGetConstraintMatrixAfterGaussElimination(CmatTmp);
        NuTo::SparseMatrixCSRVector2General<double> Cmat(CmatTmp);
        SparseMatrixCSRVector2General<double> CmatT (Cmat.Transpose());
        Eigen::VectorXd bRHS, bRHSdot, bRHSddot;
        if (CmatT.GetNumEntries() > 0)
        {
            throw MechanicsException("[NuTo::RungeKuttaBase::Solve] not implemented for constrained systems including multiple dofs.");
        }

        //calculate individual inverse mass matrix, use only lumped mass matrices - stored as fullvectors and then use asDiagonal()
        Eigen::VectorXd invMassMatrix_j(mStructure->GetNumActiveDofs());
        Eigen::VectorXd massMatrix_k;

        //extract displacements, velocities and accelerations
        Eigen::VectorXd disp_j, vel_j, tmp_k,
                                                disp_j_tmp, vel_j_tmp, //intermediate values of the displacements and velocities
                                                disp_j_new, vel_j_new; //new d and v at end of time step
        std::vector<Eigen::VectorXd > d_disp_j_tmp, d_vel_j_tmp; //intermediate values of the time derivatives of d and v

        Eigen::VectorXd extForce_j, extForce_k;
        Eigen::VectorXd outOfBalance_j(mStructure->GetNumActiveDofs()), outOfBalance_k;
        Eigen::VectorXd intForce_j(mStructure->GetNumActiveDofs()),
                                                intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

        //store last converged displacements, velocities and accelerations
        mStructure->NodeExtractDofValues(0,disp_j,tmp_k);
        mStructure->NodeExtractDofValues(1,vel_j,tmp_k);

        //calculate lumped mass matrix
        //intForce_j is just used as a tmp variable
        mStructure->BuildGlobalLumpedHession2(intForce_j,massMatrix_k);

        int contactAreas = mConstraintContact.size();
        Eigen::VectorXd tempVelAll(mStructure->GetNumActiveDofs());
        Eigen::VectorXd tempDispAll(mStructure->GetNumActiveDofs());
        if (contactAreas > 0)
        {
            for (int contact = 0; contact < contactAreas; contact++)
            {
                if (!mStructure->IsPenalty(contact))
                {
                    mStructure->SetConstraintContactMassFactor(mConstraintContact[contact], intForce_j);
                }
            }
        }

        //check the sum of all entries
        std::cout << "the total mass is " << intForce_j.sum()/((double)mStructure->GetDimension()) +  tmp_k.sum()/((double)mStructure->GetDimension()) << std::endl;

        //invert the mass matrix
        invMassMatrix_j = intForce_j.cwiseInverse();

        double curTime  = 0;
        CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
        d_disp_j_tmp.resize(this->GetNumStages());
        d_vel_j_tmp.resize(this->GetNumStages());
        std::vector<double> stageDerivativeFactor(this->GetNumStages()-1);
        while (curTime < rTimeDelta)
        {
            //calculate for delta_t = 0
            std::cout << "curTime " << curTime <<   " (" << curTime/rTimeDelta << ") max Disp = "  <<  disp_j.maxCoeff() << std::endl;
            disp_j_new = disp_j;
            vel_j_new = vel_j;

            double prevTime(mTime);
            double prevCurTime(curTime);
            for (int countStage = 0; countStage < this->GetNumStages(); countStage++)
            {
                //std::cout << "\n stage weight " << GetStageWeights(countStage) << std::endl;
                double deltaTimeStage = this->GetStageTimeFactor(countStage)*mTimeStep;
                this->GetStageDerivativeFactor(stageDerivativeFactor, countStage);
                disp_j_tmp = disp_j;
                vel_j_tmp = vel_j;

                // UNconstrained step
                for (int countStage2 = 0; countStage2 < countStage; countStage2++)
                {
                    if (stageDerivativeFactor[countStage2] != 0.)
                    {
                        disp_j_tmp+=d_disp_j_tmp[countStage2]*(stageDerivativeFactor[countStage2]);
                        vel_j_tmp +=d_vel_j_tmp [countStage2]*(stageDerivativeFactor[countStage2]);
                    }
                }

                if (this->HasTimeChanged(countStage)==true)
                {
                    curTime=prevCurTime+deltaTimeStage;
                    mTime=prevTime+deltaTimeStage;

                    //to be implemented mStructure->SetCurrentTime(mTime);
                    //an update of the external load factor and the time dependent constraint is only
                    //necessary for a modified global time
                    if (mTimeDependentConstraint!=-1)
                    {
                        throw MechanicsException("[NuTo::RungeKuttaBase::Solve] solution with constraints not yet implemented.");
                        //double timeDependentConstraintFactor = this->CalculateTimeDependentConstraintFactor(curTime);
                        //mStructure->ConstraintSetRHS(mTimeDependentConstraint,timeDependentConstraintFactor);
                        //mStructure->ConstraintGetRHSAfterGaussElimination(bRHS);
                    }
                    //calculate external force
                    this->CalculateExternalLoad(*mStructure, curTime, extForce_j, extForce_k);
                }

                // unconstrained values
                mStructure->NodeMergeActiveDofValues(0,disp_j_tmp);
                mStructure->NodeMergeActiveDofValues(1,vel_j_tmp);

                // add contact constraints
                if (contactAreas > 0 && countStage > 0)
                {
                    for (int contact = 0; contact < contactAreas; contact++)
                    {
                        if (!mStructure->IsPenalty(contact))
                        {
                            for (int countStage2 = 0; countStage2 <= countStage; countStage2++)
                            {
                                if (stageDerivativeFactor[countStage2] != 0.)
                                {
                                    disp_j_tmp+=d_disp_j_tmp[countStage2]*(stageDerivativeFactor[countStage2]);
                                    vel_j_tmp +=d_vel_j_tmp [countStage2]*(stageDerivativeFactor[countStage2]);
                                }
                            }

                            Eigen::VectorXd contactContributionVectorDisp(mStructure->GetNumActiveDofs());
                            Eigen::VectorXd contactContributionVectorVel(mStructure->GetNumActiveDofs());

                            mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVectorDisp, 0);
                            mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVectorVel, 1);

                            Eigen::VectorXd tempVel = invMassMatrix_j.asDiagonal()*contactContributionVectorVel;
                            Eigen::VectorXd tempDisp = contactContributionVectorDisp;


                            // update the rhs due to contact forces (the rhs is always multiplied by the time step)
                            d_disp_j_tmp[countStage-1]-= tempDisp/stageDerivativeFactor[countStage-1];
                            d_vel_j_tmp[countStage-1] -= tempVel/stageDerivativeFactor[countStage-1];

                            // final stage
                            disp_j_new -= tempDisp*(GetStageWeights(countStage-1))/stageDerivativeFactor[countStage-1];
                            vel_j_new  -= tempVel*(GetStageWeights(countStage-1))/stageDerivativeFactor[countStage-1];

                            // constrained step
                            disp_j_tmp -= tempDisp;
                            vel_j_tmp  -= tempVel;
                        }
                    }
                    // update displacements & velocities
                    mStructure->NodeMergeActiveDofValues(0,disp_j_tmp);
                    mStructure->NodeMergeActiveDofValues(1,vel_j_tmp);
                }

                mStructure->ElementTotalUpdateTmpStaticData();

                //calculate internal force (with update of history variables = true)
                mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j,intForce_k,false);

                //std::cout << "norm of deltaForce " << (intForce_j-extForce_j).norm() << std::endl;

                //update derivatives (ydot or k1,k2,k3,k4) for Runge Kutta
                d_disp_j_tmp[countStage] = vel_j_tmp*mTimeStep;
                //std::cout << "d_disp_j_tmp " << d_disp_j_tmp[countStage](0) << std::endl;
                d_vel_j_tmp[countStage]  = (invMassMatrix_j.asDiagonal()*(extForce_j-intForce_j))*mTimeStep;

                // Penalty
                if (contactAreas > 0 )
                {
                    for (int contact = 0; contact < contactAreas; contact++)
                    {
                        if (mStructure->IsPenalty(contact))
                        {
                            Eigen::VectorXd contactContributionVectorVel(mStructure->GetNumActiveDofs());
                            mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVectorVel, 1);

                            d_vel_j_tmp[countStage] += mTimeStep*(invMassMatrix_j.asDiagonal())*contactContributionVectorVel;
                        }
                    }
                }

                //std::cout << "d_vel_j_tmp " << d_vel_j_tmp[countStage](0) << std::endl;
                //std::cout << "norm of acc " << (d_vel_j_tmp).norm() << std::endl;

                disp_j_new += d_disp_j_tmp[countStage]*(GetStageWeights(countStage));
                //std::cout << "disp_j_new " << disp_j_new(0) << std::endl;
                vel_j_new  += d_vel_j_tmp[countStage]*(GetStageWeights(countStage));
                //std::cout << "vel_j_new " << vel_j_new(0) << std::endl;
            }

            mTime = prevTime + mTimeStep;
            curTime = prevCurTime + mTimeStep;

            //std::cout << "final disp_j_new " << disp_j_new(0) << std::endl;
            mStructure->NodeMergeActiveDofValues(0,disp_j_new);
            mStructure->NodeMergeActiveDofValues(1,vel_j_new);
            mStructure->ElementTotalUpdateTmpStaticData();
            mStructure->ElementTotalUpdateStaticData();
            //std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;
            disp_j = disp_j_new;
            vel_j  = vel_j_new;

            //**********************************************
            //PostProcessing
            //**********************************************
            if (CmatT.GetNumEntries() > 0)
            {
                throw MechanicsException("[NuTo::RungeKuttaBase::Solve] not implemented for constrained systems including multiple dofs.");
            }
            else
            {
                //outOfBalance_j is automatically zero
                //outOfBalance_j.Resize(intForce_j.rows());
                //the acceleration of the dofs k is given by the acceleration of the rhs of the constraint equation
                //this is calculated using finite differencs
                //make sure to recalculate the internal force and external force (if time factor is not 1)
                if (mTimeDependentConstraint!=-1)
                {
                    throw MechanicsException("[NuTo::RungeKuttaBase::Solve] solution with constraints not yet implemented.");
                }

                //acc_k = (bRHSprev-bRHShalf*2+bRHSend)*(4./(timeStep*timeStep))
                //outOfBalance_k = intForce_k - extForce_k + massMatrix_k.asDiagonal()*acc_k;
            }

            // postprocess data for plotting
            this->PostProcess(outOfBalance_j, outOfBalance_k);
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::RungeKuttaBase::Solve] performing Newton-Raphson iteration.");
        throw;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::HEDOPRI5Original::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::HEDOPRI5Original::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::HEDOPRI5Original::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[HEDOPRI5Original::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[HEDOPRI5Original::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[HEDOPRI5Original::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else
        {
            throw MathException ( "[Matrix::Restore]File type not implemented" );
        }
    }
    catch ( MechanicsException &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[RungeKutta4::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::RungeKutta4::Save (const std::string &filename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ofstream ofs ( filename.c_str(), std::ios_base::binary );
        std::string tmpStr ( GetTypeId() );
        if (rType=="BINARY")
        {
            boost::archive::binary_oarchive oba ( ofs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oba & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_oarchive oxa ( ofs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpStr );
            oxa & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_oarchive ota ( ofs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpStr );
            ota & boost::serialization::make_nvp(tmpStr.c_str(), *this);
        }
        else
        {
            throw MechanicsException ( "[RungeKutta4::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[RungeKutta4::Save]File save exception in boost - " ) +std::string ( e.what() ) );
        std::cout << s << "\n";
        throw MathException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw;
    }
    catch ( std::exception &e )
    {
        throw MechanicsException ( e.what() );
    }
    catch ( ... )
    {
        throw MechanicsException ( "[RungeKutta4::Save] Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::RungeKutta4)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
