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
#include "nuto/mechanics/nodes/NodeDisplacements.h"
#include "nuto/mechanics/nodes/NodeRotations.h"
#include "nuto/mechanics/nodes/NodeTemperatures.h"

// constructor
NuTo::ConstraintEquationTerm::ConstraintEquationTerm(const NodeBase* rNode, Node::eAttributes rDofType, int rDofComponent, double rCoefficient)
{
    assert(rNode != 0);
    switch (rDofType)
    {
    case Node::DISPLACEMENTS:
    {
        int numDisplacements = rNode->GetNumDisplacements();
        switch (numDisplacements)
        {
        case 1:
            if (dynamic_cast< const NuTo::NodeDisplacements<1>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of displacement dofs.");
            }
            break;
        case 2:
            if (dynamic_cast< const NuTo::NodeDisplacements<2>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of displacement dofs.");
            }
            break;
        case 3:
            if (dynamic_cast< const NuTo::NodeDisplacements<3>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of displacement dofs.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of displacement dofs.");
        }
        if (rDofComponent < 0 || rDofComponent >= numDisplacements)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::ROTATIONS:
    {
        int numRotations = rNode->GetNumRotations();
        switch (numRotations)
        {
        case 1:
            if (dynamic_cast< const NuTo::NodeRotations<1>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of rotation dofs.");
            }
            break;
        case 3:
            if (dynamic_cast< const NuTo::NodeRotations<3>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of rotation dofs.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of rotation dofs.");
        }
        if (rDofComponent < 0 || rDofComponent >= numRotations)
        {
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid displacement component.");
        }
    }
    break;
    case Node::TEMPERATURES:
    {
        int numTemperatures = rNode->GetNumTemperatures();
        switch (numTemperatures)
        {
        case 1:
            if (dynamic_cast< const NuTo::NodeTemperatures<1>* >(rNode) == 0)
            {
                throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] node does not have the correct number of temperature dofs.");
            }
            break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of temeprature dofs.");
        }
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
        int numDisplacements = this->mNode->GetNumDisplacements();
        switch (numDisplacements)
        {
        case 1:
        {
            const NuTo::NodeDisplacements<1>* tmpNodePtr = dynamic_cast< const NuTo::NodeDisplacements<1>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofDisplacement(this->mDofComponent);
        }
        break;
        case 2:
        {
            const NuTo::NodeDisplacements<2>* tmpNodePtr = dynamic_cast< const NuTo::NodeDisplacements<2>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofDisplacement(this->mDofComponent);
        }
        break;
        case 3:
        {
            const NuTo::NodeDisplacements<3>* tmpNodePtr = dynamic_cast< const NuTo::NodeDisplacements<3>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofDisplacement(this->mDofComponent);
        }
        break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of displacement dofs.");
        }
    }
    break;
    case Node::ROTATIONS:
    {
        int numRotations = this->mNode->GetNumRotations();
        switch (numRotations)
        {
        case 1:
        {
            const NuTo::NodeRotations<1>* tmpNodePtr = dynamic_cast< const NuTo::NodeRotations<1>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofRotation(this->mDofComponent);
        }
        break;
        case 3:
        {
            const NuTo::NodeRotations<3>* tmpNodePtr = dynamic_cast< const NuTo::NodeRotations<3>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofRotation(this->mDofComponent);
        }
        break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of rotation dofs.");
        }
    }
    break;
    case Node::TEMPERATURES:
    {
        int numTemperatures = this->mNode->GetNumTemperatures();
        switch (numTemperatures)
        {
        case 1:
        {
            const NuTo::NodeTemperatures<1>* tmpNodePtr = dynamic_cast< const NuTo::NodeTemperatures<1>* >(this->mNode);
            assert(tmpNodePtr != 0);
            column = tmpNodePtr->GetDofTemperature(this->mDofComponent);
        }
        break;
        default:
            throw MechanicsException("[NuTo::ConstraintEquationTerm::ConstraintEquationTerm] invalid number of temeprature dofs.");
        }
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
    ar & BOOST_SERIALIZATION_NVP(const_cast<NodeBase*&>(this->mNode))
    & BOOST_SERIALIZATION_NVP(mDofType)
    & BOOST_SERIALIZATION_NVP(mDofComponent)
    & BOOST_SERIALIZATION_NVP(mCoefficient);
}
#endif // ENABLE_SERIALIZATION
