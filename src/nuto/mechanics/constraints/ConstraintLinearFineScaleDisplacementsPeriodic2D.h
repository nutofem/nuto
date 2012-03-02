// $Id $

#ifndef CONSTRAINTFINESCALEDISPLACEMENTS2PERIODIC2D_H
#define CONSTRAINTFINESCALEDISPLACEMENTS2PERIODIC2D_H

#include "nuto/mechanics/constraints/ConstraintBase.h"
#include "nuto/mechanics/constraints/ConstraintLinear.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain2D.h"
#include "nuto/mechanics/MechanicsException.h"

namespace NuTo
{
class StructureBase;
class NodeBase;
template<class T> class Group;
class NodeCoordinatesDisplacementsMultiscale2D;
//! @author Joerg F. Unger
//! @date August 2010
//! @brief ... class for all displacement constraints periodic boundary conditions in 2D
//! in order to avoid free floating of the body, the lower left corner should be additionally fixed (to zero or any other value) in both directions
//! the general idea is a decomposition of the periodic boundary conditions into a homogeneous part (mStrain) and a discontinuous part (mCrackOpening)
//! for a crack with an angle between 45 and 235 degrees (for other angles the values are automatically adjusted by +-180
//! attention: for a tension test with a horizontal crack, either specify angle=180 and uy or angle=0 and -uy
class ConstraintLinearFineScaleDisplacementsPeriodic2D : public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLinearFineScaleDisplacementsPeriodic2D(const Group<NodeBase>* rBoundaryNodes,const EngineeringStrain2D& rStrain);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const;

    //!@brief sets/modifies the average strain applied to the boundary
    //!@param rAngle angle in deg
    void SetStrain(const EngineeringStrain2D& rStrain);

    //!@brief calculate the border vectors in counterclockwise direction
    void SetBoundaryVectors();

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void AddToConstraintMatrix(int& curConstraintEquation,
                               NuTo::SparseMatrixCSRGeneral<double>& rConstraintMatrix)const;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void GetRHS(int& curConstraintEquation,NuTo::FullMatrix<double>& rRHS)const;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short rVerboseLevel) const;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLinearFineScaleDisplacementsPeriodic2D(){};

    //! @brief calculate weighting function for first master node (linear interpolation between master nodes
    //! @param rCoordinateCurMaster
    double CalculateWeightFunction(double rCoordinateCurMaster, double rCoordinateNextMaster, double rCoordinateSlave)const;

    //! @brief calculate delta displacement in x and y direction from the applied strain and the nodal position
    void CalculateDeltaDisp(double rCoordinates[2], double rDeltaDisp[2])const;

    //! @brief average strain applied to the boundaries (epsilon_xx, epsilon_yy, gamma_xy)
    EngineeringStrain2D mStrain;

    //! @brief boundary groups (upper, lower, left and right), corner nodes are supposed to be in two groups
    const Group<NodeBase>* mGroupBoundaryNodes;

    std::vector<NodeCoordinatesDisplacementsMultiscale2D*> mSlaveNodesRightBoundary;
    std::vector<NodeCoordinatesDisplacementsMultiscale2D*> mSlaveNodesTopBoundary;
    std::vector<NodeCoordinatesDisplacementsMultiscale2D*> mMasterNodesLeftBoundary;
    std::vector<NodeCoordinatesDisplacementsMultiscale2D*> mMasterNodesBottomBoundary;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearFineScaleDisplacementsPeriodic2D)
#endif // ENABLE_SERIALIZATION

#endif //CONSTRAINTFINESCALEDISPLACEMENTS2PERIODIC2D_H
