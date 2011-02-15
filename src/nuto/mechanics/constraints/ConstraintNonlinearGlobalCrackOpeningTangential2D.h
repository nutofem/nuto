// $Id: ConstraintNode.h 314 2010-09-27 16:31:43Z unger3 $

#ifndef CONSTRAINTNONLINEARGLOBALTANGENTIALCRACKOPENING_H
#define CONSTRAINTNONLINEARGLOBALTANGENTIALCRACKOPENING_H

#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/constraints/ConstraintNonlinear.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class StructureIp;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... abstract class for all constraints applied to a single node
class ConstraintNonlinearGlobalCrackOpeningTangential2D : public ConstraintNonlinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintNonlinearGlobalCrackOpeningTangential2D(const StructureIp* rStructure, double rScalingFactor, double rPenaltyStiffness);

    ConstraintNonlinearGlobalCrackOpeningTangential2D* AsConstraintNonlinearGlobalCrackOpeningTangential2D();

    const ConstraintNonlinearGlobalCrackOpeningTangential2D* AsConstraintNonlinearGlobalCrackOpeningTangential2D()const;

    void SetPenaltyStiffness(double rPenaltyStiffness);

    void SetScalingFactor(double rScalingFactor);

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

    //! @brief calculates the coefficient matrix for the 0-th derivative in the differential equation
    //! for a mechanical problem, this corresponds to the stiffness matrix (Hessian) of the total potential
    //! @param rResult ... coefficient matrix
    //! @param rGlobalDofs ... row and column numbers in global system
    void CalculateCoefficientMatrix_0(NuTo::SparseMatrixCSRVector2Symmetric<double>& rResult,
            std::vector<int>& rGlobalDofs)const;

    //! @brief calculates the gradient of the internal potential
    //! for a mechanical problem, this corresponds to the internal force vector
    void CalculateGradientInternalPotential(NuTo::FullMatrix<double>& rResult,
            std::vector<int>& rGlobalDofs)const;

    //! @brief calculates the internal potential
    double CalculateTotalPotential()const;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintNonlinearGlobalCrackOpeningTangential2D(){};

    const StructureIp* mStructure;
    double mPenaltyStiffness;
    double mScalingFactor;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintNonlinearGlobalCrackOpeningTangential2D)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTNONLINEARGLOBALTANGENTIALCRACKOPENING_H

