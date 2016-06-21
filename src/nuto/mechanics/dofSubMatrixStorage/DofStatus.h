#pragma once



#include "nuto/mechanics/nodes/NodeEnum.h"

#ifdef ENABLE_SERIALIZATION
#include <boost/serialization/access.hpp>
#include <boost/serialization/export.hpp>
#endif // ENABLE_SERIALIZATION

#include <set>
#include <map>

namespace NuTo
{

//! @brief class to store state of the current dof setup
//! @remark gets modified by the structure after each change of the active/inactive dof types
class DofStatus
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version);
#endif // ENABLE_SERIALIZATION
public:

#ifndef SWIG

    DofStatus();

    friend std::ostream& operator<<(std::ostream& out, const DofStatus& dofStatus);

    const std::set<Node::eDof>& GetActiveDofTypes() const
    {
        return mActiveDofTypes;
    }

    void SetActiveDofTypes(const std::set<Node::eDof>& rActiveDofTypes)
    {
        mActiveDofTypes = rActiveDofTypes;
    }

    const std::set<Node::eDof>& GetDofTypes() const
    {
        return mDofTypes;
    }

    void SetDofTypes(const std::set<Node::eDof>& rDofTypes)
    {
        mDofTypes = rDofTypes;
    }

    bool HasInteractingConstraints() const
    {
        return mHasInteractingConstraints;
    }

    void SetHasInteractingConstraints(bool rHasInteractingConstraints)
    {
        mHasInteractingConstraints = rHasInteractingConstraints;
    }

    const std::map<Node::eDof, int>& GetNumActiveDofsMap() const
    {
        return mNumActiveDofs;
    }

    const std::map<Node::eDof, int>& GetNumDependentDofsMap() const
    {
        return mNumDependentDofs;
    }

    void SetNumActiveDofs(Node::eDof rDofType, unsigned int rNumActiveDofs)
    {
        mNumActiveDofs[rDofType] = rNumActiveDofs;
    }

    void SetNumDependentDofs(Node::eDof rDofType, unsigned int rNumDependentDofs)
    {
        mNumDependentDofs[rDofType] = rNumDependentDofs;
    }


    bool IsSymmetric(Node::eDof rDofType) const
    {
        return mSymmetricDofTypes.find(rDofType) != mSymmetricDofTypes.end();
    }

    void SetIsSymmetric(Node::eDof rDofType, bool rIsSymmetric)
    {
        if (rIsSymmetric)
            mSymmetricDofTypes.insert(rDofType);
        else
            mSymmetricDofTypes.erase(rDofType);
    }

#endif //SWIG

private:

    std::map<Node::eDof, int> mNumActiveDofs;
    std::map<Node::eDof, int> mNumDependentDofs;


    std::set<Node::eDof> mDofTypes;
    std::set<Node::eDof> mActiveDofTypes;

    std::set<Node::eDof> mSymmetricDofTypes;

    bool mHasInteractingConstraints;



};
} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
#ifndef SWIG
BOOST_CLASS_EXPORT_KEY(NuTo::DofStatus)
#endif // SWIG
#endif // ENABLE_SERIALIZATION
