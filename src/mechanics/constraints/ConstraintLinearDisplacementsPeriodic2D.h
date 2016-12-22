// $Id: ConstraintNodeGroupDisplacements2D.h 265 2010-06-08 08:47:00Z arnold2 $

#pragma once


#include "mechanics/constraints/ConstraintLinear.h"
#include "mechanics/MechanicsException.h"
#include "mechanics/constitutive/inputoutput/EngineeringStrain.h"

namespace NuTo
{
class StructureBase;
class NodeBase;
template<class T> class Group;
class NodeCoordinatesDisplacements2D;
//! @author Joerg F. Unger
//! @date August 2010
//! @brief ... class for all displacement constraints periodic boundary conditions in 2D
//! in order to avoid free floating of the body, the lower left corner should be additionally fixed (to zero or any other value) in both directions
//! the general idea is a decomposition of the periodic boundary conditions into a homogeneous part (mStrain) and a discontinuous part (mCrackOpening)
//! for a crack with an angle between 45 and 235 degrees (for other angles the values are automatically adjusted by +-180
//! attention: for a tension test with a horizontal crack, either specify angle=180 and uy or angle=0 and -uy
class ConstraintLinearDisplacementsPeriodic2D : public ConstraintLinear
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    //! @param rDirection ... direction of the applied constraint
    //! @param rValue ... direction of the applied constraint
    ConstraintLinearDisplacementsPeriodic2D(const StructureBase* rStructure, double rAngle, const EngineeringStrain<2>& rStrain,
            NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic> crackOpening, double rRadiusToCrackWithoutConstraints,
            const Group<NodeBase>* rGroupTop,const Group<NodeBase>* rGroupBottom,
            const Group<NodeBase>* rGroupLeft, const Group<NodeBase>* rGroupRight);

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    int GetNumLinearConstraints()const override;

    //!@brief sets/modifies angle of the boundary condition
    //!@param rAngle angle in deg
    void SetAngle(double rAngle);

    //!@brief sets/modifies the average strain applied to the boundary
    //!@param rAngle angle in deg
    void SetStrain(const EngineeringStrain<2>& rStrain) override;

    //!@brief sets/modifies the average strain applied to the boundary
    //!@param rAngle angle in deg
    void SetCrackOpening(const NuTo::FullMatrix<double,Eigen::Dynamic,Eigen::Dynamic>& rCrackOpening) override;

    //!@brief calculate the border vectors in counterclockwise direction
    void SetBoundaryVectors();

    //!@brief calculates all the nodes on the boundary
    std::vector<NodeBase*> GetBoundaryNodes();

    //! @brief adds the constraint equations to the matrix
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void AddToConstraintMatrix(int& curConstraintEquation, NuTo::SparseMatrix<double>& rConstraintMatrix)const override;

    //!@brief writes for the current constraint equation(s) the rhs into the vector
    // (in case of more than one equation per constraint, curConstraintEquation is increased based on the number of constraint equations per constraint)
    //! @param curConstraintEquation (is incremented during the function call)
    //! @param rConstraintMatrix (the first row where a constraint equation is added is given by curConstraintEquation)
    void GetRHS(int& curConstraintEquation,NuTo::FullVector<double,Eigen::Dynamic>& rRHS)const override;

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    void Info(unsigned short) const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, " to be implemented.");
    }

    //! @brief determines the dof type affected by the constraint
    //! @return dof type
    Node::eDof GetDofType() const override;

#ifdef ENABLE_SERIALIZATION
    //! @brief deserializes (loads) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void load(Archive & ar, const unsigned int version);

    //! @brief serializes (save) the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER()

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast) override;
#endif // ENABLE_SERIALIZATION

protected:
    //! @brief just for serialization
    ConstraintLinearDisplacementsPeriodic2D(){}

    //! @brief calculate weighting function for first master node (linear interpolation between master nodes
    //! @param rCoordinateCurMaster
    double CalculateWeightFunction(double rCoordinateCurMaster, double rCoordinateNextMaster, double rCoordinateSlave)const;

    //! @brief calculate delta displacement in x and y direction from the applied strain and the nodal position

    void CalculateDeltaDisp(double rCoordinates[2], double rDeltaDisp[2])const;

    //! @brief prescribed angle of the periodoc boundary conditions
    double mAngle;

    //! @brief average strain applied to the boundaries (epsilon_xx, epsilon_yy, gamma_xy)
    EngineeringStrain<2> mStrain;

    //! @brief crack opening in x and y-direction
    double mCrackOpening[2];

    //! @brief radius close to the crack, where constraints are not applied
    double mRadiusToCrackWithoutConstraints;

    //! @brief boundary groups (upper, lower, left and right), corner nodes are supposed to be in two groups
    const Group<NodeBase>* mGroupTop;
    const Group<NodeBase>* mGroupBottom;
    const Group<NodeBase>* mGroupLeft;
    const Group<NodeBase>* mGroupRight;

    const NodeBase* mLeftUpperCorner;
    const NodeBase* mLeftLowerCorner;
    const NodeBase* mRightUpperCorner;
    const NodeBase* mRightLowerCorner;

    std::vector<NodeBase*> mSlaveNodesRightBoundary;
    std::vector<NodeBase*> mSlaveNodesTopBoundary;
    std::vector<NodeBase*> mMasterNodesLeftBoundary;
    std::vector<NodeBase*> mMasterNodesBottomBoundary;

    const StructureBase* mStructure;
};
}//namespace NuTo

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintLinearDisplacementsPeriodic2D)
#endif // ENABLE_SERIALIZATION

