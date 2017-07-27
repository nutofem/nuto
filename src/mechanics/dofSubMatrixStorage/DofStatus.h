#pragma once


#include <set>
#include <map>

namespace NuTo
{

namespace Node
{
enum class eDof : unsigned char;
} // namespace Node

//! @brief class to store state of the current dof setup
//! @remark gets modified by the structure after each change of the active/inactive dof types
class DofStatus
{
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

    int GetNumActiveDofs(Node::eDof dof) const
    {
        return mNumActiveDofs.at(dof);
    }

    int GetNumDependentDofs(Node::eDof dof) const
    {
        return mNumDependentDofs.at(dof);
    }

    int GetNumDofs(Node::eDof dof) const
    {
        return GetNumActiveDofs(dof) + GetNumDependentDofs(dof);
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

    bool IsActive(Node::eDof dofType) const
    {
        return mActiveDofTypes.find(dofType) != mActiveDofTypes.end();
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

#endif // SWIG

private:
    std::map<Node::eDof, int> mNumActiveDofs;
    std::map<Node::eDof, int> mNumDependentDofs;


    std::set<Node::eDof> mDofTypes;
    std::set<Node::eDof> mActiveDofTypes;

    std::set<Node::eDof> mSymmetricDofTypes;
};
} /* namespace NuTo */
