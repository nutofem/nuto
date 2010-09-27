// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION
#include <cassert>

#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintEquationTerm.h"
#include "nuto/mechanics/nodes/NodeBase.h"

// constructor
NuTo::ConstraintEquationTerm::ConstraintEquationTerm(const NodeBase* rNode, Node::eAttributes rDofType, int rDofComponent, double rCoefficient)
{
    assert(rNode != 0);
    switch (rDofType)
    {
    case Node::DISPLACEMENTS:
    {
        int numDisplacements = rNode->GetNumDisplacements();
        if (rDofComponent < 0 || rDofComponent >= numDisplacements)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::ROTATIONS:
    {
        int numRotations = rNode->GetNumRotations();
        if (rDofComponent < 0 || rDofComponent >= numRotations)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::TEMPERATURES:
    {
        int numTemperatures = rNode->GetNumTemperatures();
        if (rDofComponent < 0 || rDofComponent >= numTemperatures)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid temerature component.");
        }
    }
    break;
    default:
        throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] dof type is not supported.");
    }

    // set values
    this->mNode = rNode;
    this->mDofType = rDofType;
    this->mDofComponent = rDofComponent;
    this->mCoefficient = rCoefficient;
}

// default constructor (should be private)
NuTo::ConstraintEquationTerm::ConstraintEquationTerm()
{
    this->mNode = 0;
    this->mDofType = Node::COORDINATES;
    this->mDofComponent = 0;
    this->mCoefficient = 0;
}

void NuTo::ConstraintEquationTerm::AddToConstraintMatrix(int rRow, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix) const
{
    assert(rRow >= 0);
    assert(rRow < rConstraintMatrix.GetNumRows());

    // determine column
    int column;
    switch (this->mDofType)
    {
    case Node::DISPLACEMENTS:
    {
        column = mNode->GetDofDisplacement(this->mDofComponent);
    }
    break;
    case Node::ROTATIONS:
    {
    	column = mNode->GetDofRotation(this->mDofComponent);
    }
    break;
    case Node::TEMPERATURES:
    {
        column = mNode->GetDofTemperature(this->mDofComponent);
    }
    break;
    default:
        throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] dof type is not supported.");
    }

    // add entry to matrix
    rConstraintMatrix.AddEntry(rRow, column, this->mCoefficient);
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquationTerm::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive> void NuTo::ConstraintEquationTerm::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintEquationTerm" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_NVP(const_cast<NodeBase*&>(this->mNode))
    & BOOST_SERIALIZATION_NVP(mDofType)
    & BOOST_SERIALIZATION_NVP(mDofComponent)
    & BOOST_SERIALIZATION_NVP(mCoefficient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintEquationTerm" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintEquationTerm)
#endif // ENABLE_SERIALIZATION
