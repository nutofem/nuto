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


#include "math/SparseMatrixCSRVector2.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "mechanics/timeIntegration/VelocityVerlet.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "mechanics/structures/Assembler.h"

#include "base/Timer.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::VelocityVerlet::VelocityVerlet (StructureBase* rStructure)  : TimeIntegrationBase (rStructure)
{
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::VelocityVerlet::Info()const
{
	TimeIntegrationBase::Info();
}

//! @brief calculate the critical time step for explicit routines
//! for implicit routines, this will simply return zero (cmp HasCriticalTimeStep())
double NuTo::VelocityVerlet::CalculateCriticalTimeStep()const
{
	double maxGlobalEigenValue = mStructure->ElementTotalCalculateLargestElementEigenvalue();
	return 2./std::sqrt(maxGlobalEigenValue);
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::VelocityVerlet::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::VelocityVerlet::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::VelocityVerlet::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::VelocityVerlet::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::VelocityVerlet::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::VelocityVerlet::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::VelocityVerlet::serialize(Archive & ar, const unsigned int version)
{
    #ifdef DEBUG_SERIALIZATION
        std::cout << "start serialization of VelocityVerlet" << "\n";
    #endif
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase);
    #ifdef DEBUG_SERIALIZATION
        std::cout << "finish serialization of VelocityVerlet" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION


//! @brief perform the time integration
//! @param rStructure ... structure
//! @param rTimeDelta ... length of the simulation
void NuTo::VelocityVerlet::Solve(double rTimeDelta)
{
#ifdef SHOW_TIME
    std::clock_t start,end;
#ifdef _OPENMP
    double wstart = omp_get_wtime ( );
#endif
    start=clock();
#endif
    NuTo::Timer timer(__FUNCTION__, mStructure->GetShowTime(), mStructure->GetLogger());

    try
   	{
   	    mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

   	    if (mStructure->GetDofStatus().HasInteractingConstraints())
   	    	throw MechanicsException(__PRETTY_FUNCTION__, "not implemented for constrained systems including multiple dofs.");

        if (mTimeStep==0.)
        {
        	if (this->HasCriticalTimeStep())
        	{
        		mTimeStep = this->CalculateCriticalTimeStep();
        	}
        	else
        	{
                throw MechanicsException("[NuTo::VelocityVerlet::Solve] time step not set for unconditional stable algorithm.");
        	}
        }

        // check constraints:
        mStructure->GetAssembler().ConstraintUpdateRhs(42);
        if (mStructure->GetAssembler().GetConstraintRhs().Export().sum() > 1.e-6)
        {
          	throw MechanicsException(__PRETTY_FUNCTION__, "solution with constraints != 0 not yet implemented.");
        }
        mStructure->GetAssembler().ConstraintUpdateRhs(0);

        std::cout << "time step " << mTimeStep << std::endl;
        std::cout << "number of time steps " << rTimeDelta/mTimeStep << std::endl;


        StructureOutputBlockVector outOfBalance(mStructure->GetDofStatus(), true);

        CalculateStaticAndTimeDependentExternalLoad();

        //store last converged displacements, velocities and accelerations
        auto dof_dt0 = mStructure->NodeExtractDofValues(0);
        auto dof_dt1 = mStructure->NodeExtractDofValues(1);
        auto dof_dt2 = mStructure->NodeExtractDofValues(2);

        //calculate lumped mass matrix
        StructureOutputBlockMatrix hessian2 = mStructure->BuildGlobalHessian2Lumped();
        double numericMass = 0;
        numericMass += hessian2.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.JK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();
        numericMass += hessian2.KK(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS).Sum();

        numericMass /= mStructure->GetDimension(); // since the mass is added to nodes in every direction
        std::cout << "the total mass is " << numericMass << std::endl;

        //invert the mass matrix
        hessian2.CwiseInvert();


        double curTime  = 0.;

        auto extLoad = CalculateCurrentExternalLoad(curTime);
        auto intForce = mStructure->BuildGlobalInternalGradient();

        while (curTime < rTimeDelta)
        {
         	//increase time step
            curTime += mTimeStep;
            if (mStructure->GetVerboseLevel() > 5)
                std::cout << "curTime " << curTime <<   " (" << curTime/rTimeDelta << ") max Disp = "  <<  dof_dt0.J[Node::eDof::DISPLACEMENTS].maxCoeff() << std::endl;

            //apply constraints for new time stepdouble RHSConstraint
//            double timeDependentConstraintFactor(0);
            //calculate new displacement approximation
            dof_dt0 += dof_dt1 * mTimeStep +  dof_dt2 * (mTimeStep*mTimeStep*0.5);

            dof_dt0.K = mStructure->NodeCalculateDependentDofValues(dof_dt0.J); //?
			mStructure->NodeMergeDofValues(0,dof_dt0);
            if (mMergeActiveDofValuesOrder1)
            {
            	dof_dt1.K = mStructure->NodeCalculateDependentDofValues(dof_dt1.J); //?
                mStructure->NodeMergeDofValues(1,dof_dt1);
            }
            if (mMergeActiveDofValuesOrder2)
            {
            	dof_dt2.K = mStructure->NodeCalculateDependentDofValues(dof_dt2.J); //?
                mStructure->NodeMergeDofValues(2,dof_dt2);
            }
			mStructure->ElementTotalUpdateTmpStaticData();
			mStructure->ElementTotalUpdateStaticData();

			//calculate external force
			extLoad = CalculateCurrentExternalLoad(curTime);

            //calculate internal force (with update of history variables = true)
			intForce = mStructure->BuildGlobalInternalGradient();

            //**********************************************
            //PostProcessing
            //**********************************************

            //postprocess data for plotting
        	this->PostProcess(extLoad-intForce);

        	//calculate new accelerations and velocities of independent dofs
        	auto dof_dt2_new = dof_dt2;
            dof_dt2_new = hessian2*(extLoad-intForce);

			dof_dt1 += (dof_dt2+dof_dt2_new) * (0.5*mTimeStep);

			//update acceleration
			dof_dt2 = dof_dt2_new;

			//update time, this is different from curTime if several load cycles are sequentially added
            mTime+=mTimeStep;
        }
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::VelocityVerlet::Solve] performing Newton-Raphson iteration.");
        throw;
    }
#ifdef SHOW_TIME
    end=clock();
#ifdef _OPENMP
    double wend = omp_get_wtime ( );
    if (mShowTime)
        mStructure->GetLogger()<<"[NuTo::VelocityVerlet::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec(" << wend-wstart <<")\n";
#else
    if (mShowTime)
        mStructure->GetLogger()<< "[NuTo::VelocityVerlet::Solve] " << difftime(end,start)/CLOCKS_PER_SEC << "sec" << "\n";
#endif
#endif
}


#ifdef ENABLE_SERIALIZATION
//! @brief ... restore the object from a file
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
//! @brief ... save the object to a file
void NuTo::VelocityVerlet::Restore (const std::string &filename, std::string rType )
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
                throw MechanicsException ( "[VelocityVerlet::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oba & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="XML")
        {
            boost::archive::xml_iarchive oxa ( ifs, std::ios::binary );
            oxa & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[VelocityVerlet::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
            oxa & boost::serialization::make_nvp(tmpString.c_str(), *this);
        }
        else if (rType=="TEXT")
        {
            boost::archive::text_iarchive ota ( ifs, std::ios::binary );
            ota & boost::serialization::make_nvp ( "Object_type", tmpString );
            if ( tmpString!=GetTypeId() )
                throw MechanicsException ( "[VelocityVerlet::Restore]Data type of object in file ("+tmpString+") is not identical to data type of object to read ("+GetTypeId() +")." );
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
        throw MechanicsException ( "[VelocityVerlet::Restore]Unhandled exception." );
    }
}

//  @brief this routine has to be implemented in the final derived classes, which are no longer abstract
//! @param filename ... filename
//! @param aType ... type of file, either BINARY, XML or TEXT
void NuTo::VelocityVerlet::Save (const std::string &filename, std::string rType )const
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
            throw MechanicsException ( "[VelocityVerlet::Save]File type not implemented." );
        }
    }
    catch ( boost::archive::archive_exception& e )
    {
        std::string s ( std::string ( "[VelocityVerlet::Save]File save exception in boost - " ) +std::string ( e.what() ) );
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
        throw MechanicsException ( "[VelocityVerlet::Save]Unhandled exception." );
    }
}

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::VelocityVerlet)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
