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

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "nuto/mechanics/constraints/ConstraintLinearNodeDisplacements1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::ConstraintLinearNodeDisplacements1D::ConstraintLinearNodeDisplacements1D(const NodeBase* rNode, double rDirection, double rValue):
        ConstraintNode(rNode), ConstraintLinear()
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintLinearNodeDisplacements1D] Length of the direction vector is zero");
    }
    // set normalized direction
    this->mDirection = rDirection / fabs(rDirection);

    // set value
    mRHS = rValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearNodeDisplacements1D::GetNumLinearConstraints()const
{
    return 1;
}

//!@brief sets/modifies the right hand side of the constraint equation
//!@param rRHS new right hand side
void NuTo::ConstraintLinearNodeDisplacements1D::SetRHS(double rRHS)
{
	mRHS=rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements1D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const
{
    // add constraint to constrain matrix
    if (mNode->GetNumDisplacements()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    }
    rConstraintMatrix.AddEntry(curConstraintEquation,mNode->GetDofDisplacement(0), this->mDirection);

    // increase constraint equation number
    curConstraintEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearNodeDisplacements1D::GetRHS(int& curConstraintEquation,NuTo::FullMatrix<double>& rRHS)const
{
    // add constraint to constrain matrix
    if (mNode->GetNumDisplacements()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintLinearNodeDisplacements1D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    }
    // set right hand side value
    rRHS(curConstraintEquation,0) = mRHS;

    // increase constraint equation number
    curConstraintEquation++;
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearNodeDisplacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLinearNodeDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearNodeDisplacements1D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRHS)
       & BOOST_SERIALIZATION_NVP(mDirection);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearNodeDisplacements1D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearNodeDisplacements1D)
BOOST_SERIALIZATION_ASSUME_ABSTRACT(NuTo::ConstraintLinearNodeDisplacements1D)
#endif // ENABLE_SERIALIZATION
