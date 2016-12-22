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

#include "math/SparseDirectSolverMUMPS.h"
#include "math/SparseDirectSolverMKLPardiso.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/groups/Group.h"
#include "mechanics/structures/StructureBase.h"
#include "mechanics/timeIntegration/NewmarkBase.h"
#include "mechanics/timeIntegration/TimeIntegrationEnum.h"
#include "math/FullMatrix.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "math/SparseMatrixCSRSymmetric.h"
#include "mechanics/structures/StructureOutputBlockVector.h"

//! @brief constructor
//! @param mDimension number of nodes
NuTo::NewmarkBase::NewmarkBase (StructureBase* rStructure)
    : TimeIntegrationBase (rStructure)
{
	mMuDampingMass = 0.;
	mToleranceForce = 1e-6;
	mMaxNumIterations = 20;
	mBeta = 0.25;
	mGamma = 0.5;
	mInternalEnergy = 0.;
	mExternalEnergy = 0.;
	mKineticEnergy = 0.;
	mDampedEnergy = 0.;

	mUseLumpedMass = false;
}

//! @brief merges the dof values depending on the numTimeDerivatives and rMergeAll
//! @param rDof_dt0 ... 0th time derivative
//! @param rDof_dt1 ... 1st time derivative
//! @param rDof_dt2 ... 2nd time derivative
//! @param rMergeAll ... false: merges dof_dt1 only when mMuMassDamping = 0, ignores dof_dt2
void NuTo::NewmarkBase::MergeDofValues(const StructureOutputBlockVector& rDof_dt0, const StructureOutputBlockVector& rDof_dt1, const StructureOutputBlockVector& rDof_dt2, bool rMergeAll)
{
    mStructure->NodeMergeDofValues(0, rDof_dt0.J, rDof_dt0.K);
    if (mStructure->GetNumTimeDerivatives() >= 1)
    {
        if (rMergeAll || mMuDampingMass == 0)
        {
            mStructure->NodeMergeDofValues(1, rDof_dt1.J, rDof_dt1.K);
        }
    }

    if (mStructure->GetNumTimeDerivatives() >= 2)
    {
        if (rMergeAll)
        {
            mStructure->NodeMergeDofValues(2, rDof_dt2.J, rDof_dt2.K);
        }
    }
    mStructure->ElementTotalUpdateTmpStaticData();
}


//! @brief ... Info routine that prints general information about the object (detail according to verbose level)
void NuTo::NewmarkBase::Info()const
{
	TimeIntegrationBase::Info();
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::NewmarkBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::NewmarkBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::NewmarkBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::NewmarkBase::serialize(Archive & ar, const unsigned int version)
{
	#ifdef DEBUG_SERIALIZATION
	    std::cout << "start serialization of Newmark" << "\n";
	#endif
	    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(TimeIntegrationBase)
	       & BOOST_SERIALIZATION_NVP(mToleranceForce)
	       & BOOST_SERIALIZATION_NVP(mMaxNumIterations)
	       & BOOST_SERIALIZATION_NVP(mBeta)
	       & BOOST_SERIALIZATION_NVP(mGamma)
	       & BOOST_SERIALIZATION_NVP(mInternalEnergy)
	       & BOOST_SERIALIZATION_NVP(mExternalEnergy)
	       & BOOST_SERIALIZATION_NVP(mKineticEnergy)
	       & BOOST_SERIALIZATION_NVP(mDampedEnergy)
	       & BOOST_SERIALIZATION_NVP(mUseLumpedMass);



    #ifdef DEBUG_SERIALIZATION
       std::cout << "finish serialization of Newmark" << "\n";
    #endif
}

#endif // ENABLE_SERIALIZATION



#ifdef ENABLE_SERIALIZATION

#ifndef SWIG
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::NewmarkBase)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
