// $Id: ConstraintLagrangeGlobalCrackAngle2D.cpp -1   $

#ifdef ENABLE_SERIALIZATION
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#endif  // ENABLE_SERIALIZATION

#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constraints/ConstraintLagrangeGlobalCrackAngle2D.h"
#include "nuto/mechanics/structures/unstructured/StructureIp.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRVector2Symmetric.h"

#define umax 1e-1

NuTo::ConstraintLagrangeGlobalCrackAngle2D::ConstraintLagrangeGlobalCrackAngle2D(const NuTo::StructureIp* rStructure) :
        ConstraintLagrange(NuTo::Constraint::EQUAL) , ConstraintBase()
{
    mStructure = rStructure;
    mLagrangeValue = 0;
    mLagrangeDOF = -1;
    mScaling = 4.*M_PI*M_PI/(umax*umax);
}

//! @brief returns the number of constraint equations
//! @return number of constraints
int NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetNumLagrangeMultipliers()const
{
    return 1;
}

//! @brief returns the Lagrange Multiplier
//! first col Lagrange, second column slack variables
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetLagrangeMultiplier(FullMatrix<double>& rLagrangeMultiplier)const
{
    rLagrangeMultiplier.Resize(1,1);
    rLagrangeMultiplier(0,0) = mLagrangeValue;
}

//! @brief sets the global dofs
//! @param rDOF current maximum DOF, this variable is increased within the routine
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::SetGlobalDofs(int& rDOF)
{
    mLagrangeDOF = rDOF++;
}

//! @brief write dof values to constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::SetGlobalDofValues(const FullMatrix<double>& rActiveDofValues, const FullMatrix<double>& rDependentDofValues)
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    double value;
    if (mLagrangeDOF >= rActiveDofValues.GetNumRows())
    {
        int dof = mLagrangeDOF - rActiveDofValues.GetNumRows();
        assert(dof < rDependentDofValues.GetNumRows());
        value = rDependentDofValues(dof,0);
    }
    else
    {
        value = rActiveDofValues(mLagrangeDOF,0);
    }
    this->mLagrangeValue = value;
}

//! @brief extract dof values from the constraints (based on global dof number)
//! @param rActiveDofValues ... active dof values
//! @param rDependentDofValues ... dependent dof values
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::GetGlobalDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    assert(rActiveDofValues.GetNumColumns() == 1);
    assert(rDependentDofValues.GetNumColumns() == 1);

    if (mLagrangeDOF >= rActiveDofValues.GetNumRows())
    {
        int dof = mLagrangeDOF - rActiveDofValues.GetNumRows();
        assert(dof < rDependentDofValues.GetNumRows());
        rDependentDofValues(dof,0) = mLagrangeValue;
    }
    else
    {
        rActiveDofValues(mLagrangeDOF,0) = mLagrangeValue;
    }
}

//! @brief renumber the global dofs according to predefined ordering
//! @param rMappingInitialToNewOrdering ... mapping from initial ordering to the new ordering
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::RenumberGlobalDofs(const std::vector<int>& rMappingInitialToNewOrdering)
{
    mLagrangeDOF = rMappingInitialToNewOrdering[mLagrangeDOF];
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackAngle2D::AsConstraintLagrange()
{
    return this;
}

//! @brief cast to Lagrange constraint - the corresponding dofs are eliminated in the global system
const NuTo::ConstraintLagrange* NuTo::ConstraintLagrangeGlobalCrackAngle2D::AsConstraintLagrange()const
{
    return this;
}

//! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
//! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
//! @param rResult ... coefficient matrix
//! @param rGlobalDofs ... row and column numbers in global system
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    //Lagrange mult + alpha+ 2 crack opening
    int dof(4);
    rResult.Resize(dof,dof);
    rGlobalDofs.resize(dof,1);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[2] = mStructure->GetDofCrackAngle();
    rGlobalDofs[3] = dofCrackOpening[0];
    rGlobalDofs[4] = dofCrackOpening[1];

    double alpha(mStructure->GetCrackAngle());
    //calculate angle orthogonal to second principal stress
    EngineeringStrain2D strainVec(mStructure->GetTotalStrain());

    FullMatrix<double> strain(2,2);
    const double *dataPtr = strainVec.GetData();
    strain(0,0) = dataPtr[0];
    strain(0,1) = dataPtr[1]*0.5;
    strain(1,0) = strain(0,1);
    strain(1,1) = dataPtr[2];

    NuTo::FullMatrix<double> eigenVectors;
    strain.EigenVectorsSymmetric(eigenVectors);
    double alpha_2 = atan2(eigenVectors(0,0),eigenVectors(1,0));

    double delta_alpha = alpha-alpha_2;
    double delta_alpha_square = delta_alpha*delta_alpha;
    double CrackOpening_square = crackOpening[0]*crackOpening[0] + crackOpening[1]*crackOpening[1];
    double g = 0.5*(delta_alpha_square-mScaling*CrackOpening_square);

    std::cout << "alpha_2" << alpha_2 << std::endl;

    if (g>-mLagrangeValue/mPenalty)
    {
        //derivative with respect to lambda and alpha
        rResult.AddEntry(0,1,delta_alpha);
        //derivative with respect to lambda and ux
        rResult.AddEntry(0,2,mScaling*crackOpening[0]);
        //derivative with respect to lambda and uy
        rResult.AddEntry(0,3,mScaling*crackOpening[1]);
        //derivative with respect to alpha^2
        rResult.AddEntry(1,1,mLagrangeValue+2.*mPenalty*(delta_alpha_square+g));
        //derivative with respect to alpha and ux
        rResult.AddEntry(1,2,-2.*mPenalty*mScaling*crackOpening[0]*delta_alpha);
        //derivative with respect to alpha and uy
        rResult.AddEntry(1,3,-2.*mPenalty*mScaling*crackOpening[1]*delta_alpha);
        //derivative with respect to ux^2
        rResult.AddEntry(2,2,-mLagrangeValue*mScaling+2.*mPenalty*mScaling*(mScaling*crackOpening[0]*crackOpening[0]-g));
        //derivative with respect to ux and uy
        rResult.AddEntry(2,3,2.*mPenalty*mScaling*mScaling*crackOpening[0]+crackOpening[1]);
        //derivative with respect to uy^2
        rResult.AddEntry(3,3,-mLagrangeValue*mScaling+2.*mPenalty*mScaling*(mScaling*crackOpening[1]*crackOpening[1]-g));
    }
    else
    {
        //derivative with respect to lambda^2
        rResult.AddEntry(0,0,-1./mPenalty);
    }
 }

//! @brief calculates the gradient of the internal potential
//! for a mechanical problem, this corresponds to the internal force vector
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
        std::vector<int>& rGlobalDofs)const
{
    int dof(4);
    rResult.Resize(dof,1);
    rGlobalDofs.resize(dof);

    boost::array<int,2> dofCrackOpening(mStructure->GetDofGlobalCrackOpening2D());
    boost::array<double,2> crackOpening(mStructure->GetGlobalCrackOpening2D());

    rGlobalDofs[0] = mLagrangeDOF;
    rGlobalDofs[2] = mStructure->GetDofCrackAngle();
    rGlobalDofs[3] = dofCrackOpening[0];
    rGlobalDofs[4] = dofCrackOpening[1];

    double alpha(mStructure->GetCrackAngle());
    //calculate angle orthogonal to second principal stress
    EngineeringStrain2D strainVec(mStructure->GetTotalStrain());

    FullMatrix<double> strain(2,2);
    const double *dataPtr = strainVec.GetData();
    strain(0,0) = dataPtr[0];
    strain(0,1) = dataPtr[1]*0.5;
    strain(1,0) = strain(0,1);
    strain(1,1) = dataPtr[2];

    NuTo::FullMatrix<double> eigenVectors;
    strain.EigenVectorsSymmetric(eigenVectors);
    double alpha_2 = atan2(eigenVectors(0,0),eigenVectors(1,0));

    double delta_alpha = alpha-alpha_2;
    double delta_alpha_square = delta_alpha*delta_alpha;
    double CrackOpening_square = crackOpening[0]*crackOpening[0] + crackOpening[1]*crackOpening[1];
    double g = 0.5*(delta_alpha_square-mScaling*CrackOpening_square);

    std::cout << "alpha_2" << alpha_2 << std::endl;

    if (g>-mLagrangeValue/mPenalty)
    {
        //derivative with respect to lambda
        rResult(0,0)= g;
        //derivative with respect to alpha
        rResult(1,0)= delta_alpha*(mLagrangeValue + 2.*mPenalty*g);
        //derivative with respect to ux
        rResult(2,0)= -mScaling*crackOpening[0]*(mLagrangeValue+2.*mPenalty*g);
        //derivative with respect to uy
        rResult(3,0)= -mScaling*crackOpening[1]*(mLagrangeValue+2.*mPenalty*g);
    }
    else
    {
        //derivative with respect to lambda
        rResult(0,0)=-mLagrangeValue/mPenalty;
    }
}

#ifdef ENABLE_SERIALIZATION
// serializes the class
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::binary_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::xml_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::text_oarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::binary_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::xml_iarchive & ar, const unsigned int version);
template void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(boost::archive::text_iarchive & ar, const unsigned int version);
template<class Archive>
void NuTo::ConstraintLagrangeGlobalCrackAngle2D::serialize(Archive & ar, const unsigned int version)
{
#ifdef DEBUG_SERIALIZATION
    std::cout << "start serialize ConstraintLagrangeGlobalCrackAngle2D" << std::endl;
#endif
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintLagrange)
       & BOOST_SERIALIZATION_BASE_OBJECT_NVP(ConstraintBase)
       & BOOST_SERIALIZATION_NVP(mStructure)
       & BOOST_SERIALIZATION_NVP(mLagrangeValue)
       & BOOST_SERIALIZATION_NVP(mLagrangeDOF)
       & BOOST_SERIALIZATION_NVP(mScaling);
#ifdef DEBUG_SERIALIZATION
    std::cout << "finish serialize ConstraintLagrangeGlobalCrackAngle2D" << std::endl;
#endif
}
BOOST_CLASS_EXPORT_IMPLEMENT(NuTo::ConstraintLagrangeGlobalCrackAngle2D)
#endif // ENABLE_SERIALIZATION
