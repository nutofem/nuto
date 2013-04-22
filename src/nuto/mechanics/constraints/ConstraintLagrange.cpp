// $Id: ConstraintLagrange.cpp 328 2010-10-01 14:39:32Z unger3 $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "nuto/math/FullVector.h"
#include "nuto/mechanics/constraints/ConstraintLagrange.h"
#include "nuto/mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstraintLagrange::ConstraintLagrange(NuTo::Constraint::eEquationSign rEquationSign) : NuTo::ConstraintNonlinear::ConstraintNonlinear()
{
    mEquationSign=rEquationSign;
    mPenalty=0.;  //standard Lagrange
}

//! @brief constructor
NuTo::ConstraintLagrange::ConstraintLagrange(NuTo::Constraint::eEquationSign rEquationSign, double rPenaltyStiffness) : NuTo::ConstraintNonlinear::ConstraintNonlinear()
{
    mEquationSign=rEquationSign;
    mPenalty=rPenaltyStiffness;
}

// destructor
NuTo::ConstraintLagrange::~ConstraintLagrange()
{
}

//! @brief sets the penalty stiffness of the augmented lagrangian
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::ConstraintLagrange::SetPenaltyStiffness(double rPenalty)
{
    mPenalty = rPenalty;
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLagrange::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrange::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrange::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrange::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrange::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrange::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLagrange::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLagrange" << std::endl;
#endif
       ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNonlinear)
          & BOOST_SERIALIZATION_NVP(mEquationSign)
          & BOOST_SERIALIZATION_NVP(mPenalty);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLagrange" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLagrange)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintLagrange)
#endif // ENABLE_SERIALIZATION
