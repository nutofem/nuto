// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/nodes/NodeDisplacements1D.h"
#include "nuto/mechanics/constraints/ConstraintNodeDisplacements1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

// constructor
NuTo::ConstraintNodeDisplacements1D::ConstraintNodeDisplacements1D(const NodeBase* rNode, double rDirection, double rValue):
        ConstraintNode(rNode)
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintNodeDisplacements1D::ConstraintNodeDisplacements1D] Length of the direction vector is zero");
    }
    // set normalized direction
    this->mDirection = rDirection / fabs(rDirection);

    // set value
    mValue = rValue;
}

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
void NuTo::ConstraintNodeDisplacements1D::SetRHS(double rRHS)
{
    mValue = rRHS;
}


//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintNodeDisplacements1D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    // set right hand side value
    rRHS(curConstraintEquation,0) = mValue;

    // add constraint to constrain matrix
    if (mNode->GetNumDisplacements()!=1)
    {
        throw MechanicsException("[NuTo::ConstraintNodeDisplacements1D::ConstraintBase] Node does not have displacements or has more than one displacement component.");
    }
    rConstraintMatrix.AddEntry(curConstraintEquation,mNode->GetDofDisplacement(0), this->mDirection);

    // increase constraint equation number
    curConstraintEquation++;
}

#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeDisplacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNodeDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNode)
    & BOOST_SERIALIZATION_NVP(mValue)
    & BOOST_SERIALIZATION_NVP(mDirection);
}
#endif // ENABLE_SERIALIZATION
