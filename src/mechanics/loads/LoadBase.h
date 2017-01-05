// $Id$
#pragma once

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <map>
#endif  // ENABLE_SERIALIZATION

#include <eigen3/Eigen/Core>
#include "mechanics/MechanicsException.h"

namespace NuTo
{
template<class T> class SparseMatrixCSRGeneral;
template <class T, int rows> class FullVector;
class NodeBase;
//! @author Jörg F. Unger, ISM
//! @date October 2009
//! @brief ... standard abstract class for all constraint equations
class LoadBase
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
#endif  // ENABLE_SERIALIZATION

public:
    //! @brief constructor
    LoadBase(){}

    //! @brief constructor
    LoadBase(int rLoadCase);

    //! @brief ... destructor
    virtual ~LoadBase(){}

    //! @brief adds the load to global sub-vectors
    //! @param rLoadCase number of the current load case
    //! @param rActiceDofsLoadVector ... global load vector which correspond to the active dofs
    //! @param rDependentDofsLoadVector ... global load vector which correspond to the dependent dofs
    virtual void AddLoadToGlobalSubVectors(int rLoadCase, Eigen::VectorXd& rActiceDofsLoadVector, Eigen::VectorXd& rDependentDofsLoadVector)const=0;


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(mLoadCase);
    }

    //! @brief NodeBase-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetNodePtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mNodeMapCast)
    {
        (void)mNodeMapCast;
        /* Do nothing until needed, see e.g. LoadNode-class*/
    }


    //! @brief Element-Pointer are not serialized to avoid cyclic dependencies, but are serialized as Pointer-Address (uintptr_t)
    //! Deserialization of the NodeBase-Pointer is done by searching and casting back the Address in the map
    //! @param mNodeMapCast   std::map containing the old and new Addresses
    virtual void SetElementPtrAfterSerialization(const std::map<std::uintptr_t, std::uintptr_t>& mElementMapCast)
    {
        (void)mElementMapCast;
        /* Do nothing until needed, see e.g. LoadSurfaceBase2D-class*/
    }
#endif // ENABLE_SERIALIZATION


protected:
    int mLoadCase;

};
}//namespace NuTo
#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::LoadBase)
#endif // ENABLE_SERIALIZATION

