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
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/constraints/ConstraintNodeGroupDisplacements1D.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

//! @brief constructor
NuTo::ConstraintNodeGroupDisplacements1D::ConstraintNodeGroupDisplacements1D(const Group<NodeBase>* rGroup, double rDirection, double rValue) :
        ConstraintNodeGroup(rGroup)
{
    // set direction
    if (fabs(rDirection) < 1e-14)
    {
        throw MechanicsException("[NuTo::ConstraintNodeGroupDisplacements1D::ConstraintNodeGroupDisplacements1D] Length of the direction vector is zero");
    }
    // set normalized direction
    this->mDirection = rDirection / fabs(rDirection);

    // set value
    mValue = rValue;
}

//!@brief sets/modifies the right hand side of the constraint equations
//!@param rRHS new right hand side
void NuTo::ConstraintNodeGroupDisplacements1D::SetRHS(double rRHS)
{
    mValue = rRHS;
}

//! @brief adds the constraint equations to the matrix
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
//! @param rRHS right hand side of the constraint equation
void NuTo::ConstraintNodeGroupDisplacements1D::AddToConstraintMatrix(int& curConstraintEquation,
        NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,
        NuTo::FullMatrix<double>& rRHS)const
{
    // loop over nodes
    for (Group<NodeBase>::const_iterator itNode=mGroup->begin(); itNode!=mGroup->end(); itNode++)
    {
        // set right hand side value
        rRHS(curConstraintEquation,0) = mValue;

        // add constraint to constrain matrix
        if ((*itNode)->GetNumDisplacements()!=1)
        {
            throw MechanicsException("[NuTo::ConstraintNodeGroupDisplacements1D::AddToConstraintMatrix] Node does not have displacements or has more than one displacement component.");
        }
        rConstraintMatrix.AddEntry(curConstraintEquation,(*itNode)->GetDofDisplacement(0),this->mDirection);

        // increase constraint equation number
        curConstraintEquation++;
    }
}


#ifdef ENABLE_SERIALIZATION
// serialize
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintNodeGroupDisplacements1D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintNodeGroupDisplacements1D::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintNodeGroup)
    & BOOST_SERIALIZATION_NVP(mValue)
    & BOOST_SERIALIZATION_NVP(mDirection);
}
#endif // ENABLE_SERIALIZATION
