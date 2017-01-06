#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#endif  // ENABLE_SERIALIZATION

#include <iostream>

#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constraints/ConstraintEquationTerm.h"
#include "mechanics/constraints/ConstraintLinearEquation.h"
#include "mechanics/nodes/NodeBase.h"
#include "mechanics/nodes/NodeEnum.h"

// constructor
NuTo::ConstraintLinearEquation::ConstraintLinearEquation(const NodeBase* rNode, Node::eDof rDofType, int rDofComponent, double rCoefficient, double rRhsValue)
{
    this->AddTerm(rNode, rDofType, rDofComponent, rCoefficient);
    this->mRhsValue = rRhsValue;
}

// add term
void NuTo::ConstraintLinearEquation::AddTerm(const NodeBase* rNode, Node::eDof rDofType, int rDofComponent, double rCoefficient)
{
    try
    {
        if (mTerms.size() != 0)
            if (mTerms[0].GetDofType() != rDofType)
            {
                throw NuTo::MechanicsException(std::string("[NuTo::ConstraintLinearEquation::AddTerm] Constraints between different dof types (")
                + Node::DofToString(mTerms[0].GetDofType()) + " and " + Node::DofToString(rDofType) + ") currently not supported. Please (pleeeaaaaase!!!), feel free to implement it.");
            }

        ConstraintEquationTerm term(rNode, rDofType, rDofComponent, rCoefficient);
        this->mTerms.push_back(term);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::ConstraintLinearEquation::AddTerm] error creating a new constraint term");
        throw;
    }
}

// add constraint equation to constraint matrix
void NuTo::ConstraintLinearEquation::AddToConstraintMatrix(int& rConstraintLinearEquation,
		NuTo::SparseMatrix<double>& rConstraintMatrix) const
{
    // loop over constraint terms
    for (unsigned int termCount = 0; termCount < this->mTerms.size(); termCount++)
    {
        // add terms to constraint matrix
        this->mTerms[termCount].AddToConstraintMatrix(rConstraintLinearEquation, rConstraintMatrix);
    }

    // increase constraint equation number
    rConstraintLinearEquation++;
}

//!@brief writes for the current constraint equation(s) the rhs into the vector
// (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
//! @param curConstraintEquation (is incremented during the function call)
//! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
void NuTo::ConstraintLinearEquation::GetRHS(int& rConstraintLinearEquation, Eigen::VectorXd& rRHS)const
{
    // set right hand side value
    rRHS(rConstraintLinearEquation) = this->mRhsValue;

    // increase constraint equation number
    rConstraintLinearEquation++;
}

//!@brief returns the rhs
double NuTo::ConstraintLinearEquation::GetRHS()const
{
	return this->mRhsValue;
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLinearEquation::GetNumLinearConstraints()const
{
    return 1;
}

NuTo::Node::eDof NuTo::ConstraintLinearEquation::GetDofType() const
{
    if (mTerms.size() != 0)
        return mTerms[0].GetDofType();
    else
        throw MechanicsException("[NuTo::ConstraintLinearEquation::GetDofType] Equation has no terms.");
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLinearEquation::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive> void NuTo::ConstraintLinearEquation::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLinearEquation" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLinear)
       & BOOST_SERIALIZATION_NVP(mRhsValue)
       & BOOST_SERIALIZATION_NVP(mTerms);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLinearEquation" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLinearEquation)
#endif // ENABLE_SERIALIZATION
