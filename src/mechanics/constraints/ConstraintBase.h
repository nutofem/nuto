// $Id$

#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <map>
#endif  // ENABLE_SERIALIZATION


#include <eigen3/Eigen/Dense>



namespace NuTo
{
class NodeBase;
class ConstraintLinear;
class ConstraintNonlinear;
class ConstraintLagrange;
template<int TDim> class EngineeringStrain;
template<class T> class SparseMatrixCSRGeneral;




namespace Node
{
    enum class eDof : unsigned char;
}// namespace Node

//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations
class ConstraintBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    ConstraintBase();

    //! @brief destructor
    virtual ~ConstraintBase();

    //! @brief true, if the constraint is linear
    virtual bool IsLinear()const=0;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumLinearConstraints()const;

    //! @brief returns the number of constraint equations
    //! @return number of constraints
    virtual int GetNumLagrangeMultipliers()const;

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    virtual ConstraintLinear* AsConstraintLinear();

    //! @brief cast to linear constraint - the corresponding dofs are eliminated in the global system
    virtual const ConstraintLinear* AsConstraintLinear()const;

    //! @brief cast to nonlinear constraint - the corresponding dofs are eliminated in the global system
    virtual ConstraintNonlinear* AsConstraintNonlinear();

    //! @brief cast to nonlinear constraint - the corresponding dofs are eliminated in the global system
    virtual const ConstraintNonlinear* AsConstraintNonlinear()const;

    //!@brief sets/modifies the right hand side of the constraint equations
    //!@param rRHS new right hand side
    virtual void SetRHS(double rRHS);

    //!@brief returns the right hand side of the constraint equations
    //!@return rRHS
    virtual double GetRHS()const;

    //!@brief set the strain of the periodic boundary conditions
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual void SetStrain(const EngineeringStrain<2>& rStrain);

    //!@brief get the strain of a constrain equation
    //!@param rStrain strain (e_xx,e_yy,gamma_xy)
    virtual const EngineeringStrain<2>& GetStrain()const;

    //! @brief exchanges the node ptr in the full data set (elements, groups, loads, constraints etc.)
    //! this routine is used, if e.g. the data type of a node has changed, but the restraints, elements etc. are still identical
    virtual void ExchangeNodePtr(NodeBase* rOldPtr, NodeBase* rNewPtr){}

    //! @brief ... print information about the object
    //! @param rVerboseLevel ... verbosity of the information
    virtual void Info(unsigned short rVerboseLevel) const=0;

    //! @brief determines the dof type affected by the constraint
    //! @return dof type
    virtual Node::eDof GetDofType() const = 0;

#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
    {
        (void)mNodeMapCast;
        /* Do nothing until needed, see e.g. ConstraintNode-class*/
    }

    //! @brief ElementBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the ElementBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast std::map containing the old and new Addresses
    virtual void SetElementPtrAfterSerialization(const std::map<uintptr_t, uintptr_t>& mElementMapCast)
    {
        (void)mElementMapCast;
        /* Do nothing until needed, see e.g. ConstraintLinearDerivativeNonlocalTotalStrain1D-class*/
    }

#endif // ENABLE_SERIALIZATION


protected:
};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::ConstraintBase)
#endif // ENABLE_SERIALIZATION
