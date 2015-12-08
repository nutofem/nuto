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

#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/timeIntegration/ExplicitEuler.h"
#include "nuto/mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "nuto/math/FullMatrix.h"


NuTo::ExplicitEuler::ExplicitEuler (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
    mTimeStep = 0.;
    mOverallTime = 0.;
    mTimeComputation = 0.;
}

//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::ExplicitEuler::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
//! this is the critical time step from velocity verlet, the real one is certainly larger
double NuTo::ExplicitEuler::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
    return 2.0/std::sqrt(maxGlobalEigenValue);
}


#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ExplicitEuler::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ExplicitEuler::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ExplicitEuler::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ExplicitEuler::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ExplicitEuler::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ExplicitEuler::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ExplicitEuler::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of ExplicitEuler" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ExplicitEuler);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of ExplicitEuler" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief ... Return the name of the class, this is important for the serialize routines, since this is stored in the ExplicitEulere
//!            in case of restoring from a ExplicitEulere with the wrong object type, the ExplicitEulere id is printed
//! @return    class name
std::string NuTo::ExplicitEuler::GetTypeId()const
{
    return std::string("ExplicitEuler");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of ExplicitEuler" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
           & BOOST_SERIALIZATION_NVP(mTimeStep);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of ExplicitEuler" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
NuTo::Error::eError NuTo::ExplicitEuler::Solve(double rTimeDelta)
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
                throw MechanicsException("[NuTo::ExplicitEuler::Solve] time step not set for unconditional stable algorithm.");
            }
        }

        std::cout << "modify computation of critical time step to include the dependence on the time integration scheme." << std::endl;
        //calculate instead the smallest eigenfrequency, depending on the time integration this gives the critical time step

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;

        //renumber dofs and build constraint matrix
        mStructure->NodeBuildGlobalDofs();

        //calculate individual inverse mass matrix, use only lumped mass matrices - stored as fullvectors and then use asDiagonal()
        NuTo::FullVector<double,Eigen::Dynamic> invMassMatrix_j(mStructure->GetNumActiveDofs());
        NuTo::FullVector<double,Eigen::Dynamic> massMatrix_k;

        //extract displacements, velocities and accelerations
        NuTo::FullVector<double,Eigen::Dynamic> disp_j, vel_j, tmp_k,
                                                disp_j_new, vel_j_new; //new d and v at end of time step

        NuTo::FullVector<double,Eigen::Dynamic> extForce_j, extForce_k;
        NuTo::FullVector<double,Eigen::Dynamic> outOfBalance_j(mStructure->GetNumActiveDofs()), outOfBalance_k;
        NuTo::FullVector<double,Eigen::Dynamic> intForce_j(mStructure->GetNumActiveDofs()),
                                                intForce_k(mStructure->GetNumDofs() - mStructure->GetNumActiveDofs());

        //store last converged displacements, velocities and accelerations
        mStructure->NodeExtractDofValues(0, disp_j, tmp_k);
        mStructure->NodeExtractDofValues(1, vel_j , tmp_k);

        //calculate lumped mass matrix
        //intForce_j is just used as a tmp variable
        mStructure->BuildGlobalLumpedHession2(intForce_j,massMatrix_k);

        int contactAreas = mConstraintContact.size();
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

        mTimeComputation = 0.;
        mOverallTime = 0.;

#ifdef _OPENMP
    double startWtimeAll = omp_get_wtime();
#endif
        while (curTime < rTimeDelta)
        {
#ifdef _OPENMP
    double startWtime = omp_get_wtime();
#endif
            //calculate for delta_t = 0
            std::cout << "curTime " << curTime <<   " (" << curTime/rTimeDelta << ") max Disp = "  <<  disp_j.maxCoeff() << std::endl;

            double prevTime(mTime);
            double prevCurTime(curTime);

            ///////////////////// calculate forces /////////////////////
            // calculate internal force (with update of history variables = true)
            mStructure->BuildGlobalGradientInternalPotentialSubVectors(intForce_j, intForce_k, true);

            /////////////////////// calculate displacements and velocities /////////////////////////
            // unconstrained step
            disp_j_new = disp_j + mTimeStep*vel_j;
            vel_j_new = vel_j - mTimeStep*invMassMatrix_j.asDiagonal()*intForce_j;

            mStructure->NodeMergeActiveDofValues(0,disp_j_new);

            // lagrange constarained step
            for (int contact = 0; contact < contactAreas; contact++)
            {
                if (!mStructure->IsPenalty(contact))
                {
                    NuTo::FullVector<double,Eigen::Dynamic> contactContributionVectorDisp(mStructure->GetNumActiveDofs());
                    NuTo::FullVector<double,Eigen::Dynamic> contactContributionVectorVel(mStructure->GetNumActiveDofs());

                    mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVectorDisp, 0);
                    mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVectorVel, 1);

                    disp_j_new -= contactContributionVectorDisp;
                    vel_j_new  -= (invMassMatrix_j.asDiagonal())*contactContributionVectorVel;
                }
            }

            // Penalty
            if (contactAreas > 0 )
            {
                for (int contact = 0; contact < contactAreas; contact++)
                {
                    if (mStructure->IsPenalty(contact))
                    {
                        NuTo::FullVector<double,Eigen::Dynamic> contactContributionVector(mStructure->GetNumActiveDofs());
                        mStructure->SetConstraintContactVector(mConstraintContact[contact], contactContributionVector, 1);

                        vel_j_new += mTimeStep*(invMassMatrix_j.asDiagonal())*contactContributionVector;
                    }
                }
            }

            mTime = prevTime + mTimeStep;
            curTime = prevCurTime + mTimeStep;

            mStructure->NodeMergeActiveDofValues(0,disp_j_new);
            mStructure->NodeMergeActiveDofValues(0,vel_j_new);
            mStructure->ElementTotalUpdateTmpStaticData();
            mStructure->ElementTotalUpdateStaticData();
            //std::cout << "delta disp between time steps" <<  (disp_j-disp_j_new).norm() << std::endl;
            disp_j = disp_j_new;
            vel_j = vel_j_new;

#ifdef _OPENMP
    double endWtime = omp_get_wtime();
    mTimeComputation += endWtime - startWtime;
#endif
            // postprocess data for plotting
            this->PostProcess(outOfBalance_j, outOfBalance_k);
        }

#ifdef _OPENMP
    double endWtimeAll = omp_get_wtime();
    mOverallTime = endWtimeAll - startWtimeAll;
#endif

    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::ExplicitEuler::Solve] performing Newton-Raphson iteration.");
        throw e;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::ExplicitEuler::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::ExplicitEuler::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
    return NuTo::Error::SUCCESSFUL;

}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a ExplicitEulere
//! @param ExplicitEulerename ... ExplicitEulerename
//! @param aType ... type of ExplicitEulere, either BINARY, XML or TEXT
//! @brief ... save the object to a ExplicitEulere
void NuTo::ExplicitEuler::Restore (const std::string &ExplicitEulerename, std::string rType )
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ifstream ifs ( ExplicitEulerename.c_str(), std::ios_base::binary );
        std::string tmpString;
        if (rType=="BINARY")
        {
            boost::archive::binary_iarchive oba ( ifs, std::ios::binary );
            oba & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[ExplicitEuler::Restore]Data type of object in ExplicitEulere ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[ExplicitEuler::Restore]Data type of object in ExplicitEulere ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[ExplicitEuler::Restore]Data type of object in ExplicitEulere ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            ota & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else
        {
            throw MathException ( "[Matrix::Restore]ExplicitEulere type not implemented" );
        }
    }
    catch ( MechanicsException &e )
    {
        throw e;
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
//! @param ExplicitEulerename ... ExplicitEulerename
//! @param aType ... type of ExplicitEulere, either BINARY, XML or TEXT
void NuTo::RungeKutta4::Save (const std::string &ExplicitEulerename, std::string rType )const
{
    try
    {
        //transform to uppercase
        std::transform(rType.begin(), rType.end(), rType.begin(), toupper);
        std::ofstream ofs ( ExplicitEulerename.c_str(), std::ios_base::binary );
        std::string tmpStr ( GetTypeId() );
        std::string baseClassStr = tmpStr.substr ( 4,100 );
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
            throw MechanicsException ( "[RungeKutta4::Save]ExplicitEulere type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception e )
    {
        std::string s ( std::string ( "[RungeKutta4::Save]ExplicitEulere save exception in boost - " ) +std::string ( e.what() ) );
        std::cout << s << "\n";
        throw MathException ( s );
    }
    catch ( MechanicsException &e )
    {
        throw e;
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
