// $Id$

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintEquation.h"

// constructor
NuTo::ConstraintEquation::ConstraintEquation(const NodeBase* rNode, NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient, double rRhsValue)
{
    this->AddTerm(rNode, rDofType, rDofComponent, rCoefficient);
    this->mRhsValue = rRhsValue;
}

// add term
void NuTo::ConstraintEquation::AddTerm(const NodeBase* rNode, NodeBase::eAttributes rDofType, int rDofComponent, double rCoefficient)
{
    try
    {
        ConstraintEquationTerm term(rNode, rDofType, rDofComponent, rCoefficient);
        this->mTerms.push_back(term);
    }
    catch (MechanicsException& e)
    {
        e.AddMessage("[NuTo::ConstraintEquation::AddTerm] error creating a new constraint term");
        throw e;
    }
}

// add constraint equation to constraint matrix
void NuTo::ConstraintEquation::AddToConstraintMatrix(int& rConstraintEquation, NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix,NuTo::FullMatrix<double>& rRHS) const
{
    // loop over constraint terms
    for (unsigned int termCount = 0; termCount < this->mTerms.size(); termCount++)
    {
        // add terms to constraint matrix
        this->mTerms[termCount].AddToConstraintMatrix(rConstraintEquation, rConstraintMatrix);
    }

    // set right hand side value
    rRHS(rConstraintEquation,0) = this->mRhsValue;

    // increase constraint equation number
    rConstraintEquation++;
}

// get num constraint equations
int NuTo::ConstraintEquation::GetNumConstraintEquations() const
{
    return 1;
}

// set right-hand-side
void NuTo::ConstraintEquation::SetRHS(double rRHS)
{
    this->mRhsValue = rRHS;
}

#ifdef ENABLE_SERIALIZATION
template void NuTo::ConstraintEquation::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquation::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquation::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquation::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquation::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintEquation::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive> void NuTo::ConstraintEquation::serialize(Archive & ar, const unsigned int version)
{
    ar & BOOST_SERIALIZATION_NVP(mTerms)
    & BOOST_SERIALIZATION_NVP(mRhsValue);
}
#endif // ENABLE_SERIALIZATION
