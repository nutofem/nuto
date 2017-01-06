// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "mechanics/constraints/ConstraintBase.h"
#include "mechanics/constraints/ConstraintLagrange.h"
#include "mechanics/constraints/ConstraintLinear.h"
#include "mechanics/MechanicsException.h"

//! @brief constructor
NuTo::ConstraintBase::ConstraintBase()
{
}

// destructor
NuTo::ConstraintBase::~ConstraintBase()
{
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintBase::GetNumLinearConstraints()const
{
    return 0;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintBase::GetNumLagrangeMultipliers()const
{
    return 0;
}

//! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintLinear* NuTo::ConstraintBase::AsConstraintLinear()
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint is not linear.");
}

//! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintLinear* NuTo::ConstraintBase::AsConstraintLinear()const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint is not linear.");
}

//! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintNonlinear* NuTo::ConstraintBase::AsConstraintNonlinear()
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint is not nonlinear.");
}

//! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintNonlinear* NuTo::ConstraintBase::AsConstraintNonlinear()const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint is not nonlinear.");
}

//! @brief cast to linear constraint - Lagrange multipliers are added to the system of equations
NuTo::ConstraintLagrange* NuTo::ConstraintBase::AsConstraintLagrange()
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint has no Lagrange multipliers.");
}

//! @brief cast to linear constraint - Lagrange multipliers are added to the system of equations
const NuTo::ConstraintLagrange* NuTo::ConstraintBase::AsConstraintLagrange()const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Constraint has no Lagrange multipliers.");
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintBase::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintBase::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintBase::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintBase" << std::endl;
#endif
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintBase" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintBase)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintBase)
#endif // ENABLE_SERIALIZATION


//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
void NuTo::ConstraintBase::SetRHS(double rRHS)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Set right hand side for this type of constraints not implemented.");
}

//!@brief returns the right hand side of the constraint equations
//!@return rRHS
double NuTo::ConstraintBase::GetRHS()const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Get right hand side for this type of constraints not implemented.");
}

//!@brief set the strain of the periodic boundary conditions
//!@param rStrain strain (e_xx,e_yy,gamma_xy)
void NuTo::ConstraintBase::SetStrain(const NuTo::EngineeringStrain<2>& rStrain)
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Set strain for this type of constraints not implemented.");
}

//!@brief get the strain of the strain for a constrain equation
//!@return rStrain strain (e_xx,e_yy,gamma_xy)
const NuTo::EngineeringStrain<2>& NuTo::ConstraintBase::GetStrain()const
{
    throw MechanicsException(__PRETTY_FUNCTION__, "Get strain for this type of constraints not implemented.");
}


