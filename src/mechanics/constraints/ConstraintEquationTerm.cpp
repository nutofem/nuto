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
#include <iostream>


#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constraints/ConstraintEquationTerm.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

// constructor
NuTo::ConstraintEquationTerm::ConstraintEquationTerm(const NodeBase* rNode, Node::eDof rDofType, int rDofComponent, double rCoefficient)
{
    assert(rNode != 0);
    switch (rDofType)
    {
    case Node::eDof::DISPLACEMENTS:
    {
        int numDisplacements = rNode->GetNum(Node::eDof::DISPLACEMENTS);
        if (rDofComponent < 0 || rDofComponent >= numDisplacements)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::eDof::ROTATIONS:
    {
        int numRotations = rNode->GetNum(Node::eDof::ROTATIONS);
        if (rDofComponent < 0 || rDofComponent >= numRotations)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::eDof::TEMPERATURE:
    {
        int numTemperature = rNode->GetNum(Node::eDof::TEMPERATURE);
        if (rDofComponent < 0 || rDofComponent >= numTemperature)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid temerature component.");
        }
    }
    break;
    case Node::eDof::NONLOCALEQSTRAIN:
    {
        int numNonlocalEqStrain = rNode->GetNum(Node::eDof::NONLOCALEQSTRAIN);
        if (rDofComponent < 0 || rDofComponent >= numNonlocalEqStrain)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid nonlocal eq strain component.");
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

void NuTo::ConstraintEquationTerm::AddToConstraintMatrix(int rRow, NuTo::SparseMatrix<double>& rConstraintMatrix) const
{
    assert(rRow >= 0);
    assert(rRow < rConstraintMatrix.GetNumRows());

    // add entry to matrix
    rConstraintMatrix.AddValue(rRow, this->GetDof(), this->mCoefficient);
}

int NuTo::ConstraintEquationTerm::GetDof() const
{
    // determine dof
    switch (this->mDofType)
    {
    case Node::eDof::DISPLACEMENTS:
    {
        return mNode->GetDof(Node::eDof::DISPLACEMENTS, this->mDofComponent);
    }
    break;
    case Node::eDof::ROTATIONS:
    {
        return mNode->GetDof(Node::eDof::ROTATIONS, this->mDofComponent);
    }
    break;
    case Node::eDof::TEMPERATURE:
    {
        return mNode->GetDof(Node::eDof::TEMPERATURE, 0);
    }
    break;
    case Node::eDof::NONLOCALEQSTRAIN:
    {
        return mNode->GetDof(Node::eDof::NONLOCALEQSTRAIN, 0);
    }
    break;
    default:
        throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] dof type is not supported.");
    }
}

#ifdef ENABLE_SERIALIZATION

NuTo::ConstraintEquationTerm::ConstraintEquationTerm()
{
    this->mNode = 0;
    this->mDofType = Node::eDof::COORDINATES;
    this->mDofComponent = 0;
    this->mCoefficient = 0;
}

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
    ar & BOOST_SERIALIZATION_NVP(mNode)
       & BOOST_SERIALIZATION_NVP(mDofType)
       & BOOST_SERIALIZATION_NVP(mDofComponent)
       & BOOST_SERIALIZATION_NVP(mCoefficient);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintEquationTerm" << std::endl;
#endif
}

BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintEquationTerm)
#endif // ENABLE_SERIALIZATION
