#include "mechanics/dofSubMatrixStorage/BlockScalar.h"
#include "mechanics/dofSubMatrixStorage/DofStatus.h"

#include "mechanics/nodes/NodeEnum.h"
#include <cassert>
#include <iostream>
#include <iomanip>

//! @brief constructor
//! @param rDofStatus ... reference to DofStatus
NuTo::BlockScalar::BlockScalar(const NuTo::DofStatus &rDofStatus)
    : BlockStorageBase(rDofStatus)
{
}




//! @brief copy assignment operator
//! @param rOther ... other BlockScalar
NuTo::BlockScalar& NuTo::BlockScalar::operator=(const NuTo::BlockScalar &rOther)
{
    mData = rOther.mData;
    return *this;
}



//! @brief == operator
//! @param rOther ... other BlockScalar
bool NuTo::BlockScalar::operator==(const NuTo::BlockScalar &rOther) const
{
    for(auto dof : mDofStatus.GetDofTypes())
    {
        if((*this)[dof]!= rOther[dof])
        {
            return false;
        }
    }
    return true;
}

//! @brief != operator
//! @param rOther ... other BlockScalar
bool NuTo::BlockScalar::operator!=(const NuTo::BlockScalar &rOther) const
{
    return !((*this)==rOther);
}



//! @brief [] operator
//! @param rDof ... degree of freedom
//! @return reference to map entry, adds the entry to map, if needed
double &NuTo::BlockScalar::operator[](Node::eDof rDof)
{
    return mData[rDof];
}

//! @brief [] operator
//! @param rDof ... degree of freedom
//! @return reference to map entry
const double NuTo::BlockScalar::operator[](NuTo::Node::eDof rDof) const
{
    assert(mData.find(rDof)!=mData.end());
    return mData.at(rDof);
}



//! @brief *= operator
//! @param rRhs ... right operand
NuTo::BlockScalar& NuTo::BlockScalar::operator*=(double rRhs)
{
    for(auto dof : mDofStatus.GetDofTypes())
    {
        mData[dof]*=rRhs;
    }
    return *this;
}


//! @brief /= operator (only active DOFs)
//! @param rRhs ... right operand
NuTo::BlockScalar& NuTo::BlockScalar::operator/=(double rRhs)
{
    for(auto dof : mDofStatus.GetDofTypes())
    {
        mData[dof]/=rRhs;
    }
    return *this;
}

namespace NuTo
{

std::ostream& operator<<(std::ostream& rOstream, const BlockScalar& rRhs)
{
    std::ios::fmtflags f(rOstream.flags());
    rOstream << std::scientific << std::setprecision(3);
    constexpr int dofNameLength = 4;
    for (auto dof : rRhs.GetDofStatus().GetActiveDofTypes())
    {
        std::string dofName = Node::DofToString(dof);
        if (dofName.length() > dofNameLength)
            dofName = dofName.substr(0, dofNameLength);
        rOstream << "[" << dofName << " : " << rRhs[dof] << "] ";
    }
    rOstream.flags(f);
    return rOstream;
}

}  // namespace NuTo


//! @brief defines a default value to all uninitialized dof types
//! @param rDefaultValue ... default value
void NuTo::BlockScalar::DefineDefaultValueToIninitializedDofTypes(double rDefaultValue)
{
    for (auto dof : mDofStatus.GetDofTypes())
        if (mData.find(dof) == mData.end())
            mData[dof] = rDefaultValue;
}

//! @brief Check if each active DOF value of this BlockScalar is smaller than the DOF value from the other BlockScalar
//! @param rOther ... other BlockScalar
//! @return true or false
bool NuTo::BlockScalar::CheckDofWiseLessActivDofs(const NuTo::BlockScalar &rOther) const
{

    for(auto dof : mDofStatus.GetActiveDofTypes())
    {
        if((*this)[dof]>=rOther[dof])
        {
            return false;
        }
    }
    return true;
}


//! @brief Check if each active DOF value of this BlockScalar is greater than the DOF value from the other BlockScalar
//! @param rOther ... other BlockScalar
//! @return true or false
bool NuTo::BlockScalar::CheckDofWiseGreaterActivDofs(const NuTo::BlockScalar &rOther) const
{
    for(auto dof : mDofStatus.GetActiveDofTypes())
    {
        if((*this)[dof]<=rOther[dof])
        {
            return false;
        }
    }
    return true;
}



//! @brief gets the number of columns of the block storage for a specific set of dofs
//! @param rDofTypes ... set of dofs
//! @return number of columns
int NuTo::BlockScalar::GetNumColumnsDof(const std::set<NuTo::Node::eDof> &rDofTypes) const
{
    return 1;
}

//! @brief gets the number of rows of the block storage for a specific set of dofs
//! @param rDofTypes ... set of dofs
//! @return number of rows
int NuTo::BlockScalar::GetNumRowsDof(const std::set<NuTo::Node::eDof> &rDofTypes) const
{
#ifdef DEBUG
    for(auto dof : rDofTypes)
    {
        assert(mData.find(dof)!=mData.end());
    }
#endif
    return rDofTypes.size();
}


//! @brief info method
void NuTo::BlockScalar::Info() const
{
    std::cout << std::endl;
    std::cout << "BlockScalar values" << std::endl;
    std::cout << "------------------" << std::endl;
    for(auto it_Values : mData)
    {
        std::string dofName = Node::DofToString(it_Values.first);
        std::string dofStatus = "(  active  )";
        int numAddtionalBlanks = 25-dofName.length();
        if(mDofStatus.GetActiveDofTypes().find(it_Values.first) == mDofStatus.GetActiveDofTypes().end())
        {
            dofStatus = "( inactive )";
        }
        std::cout << dofName << std::string(numAddtionalBlanks, '.') << dofStatus << ": " << it_Values.second << std::endl;
    }
    std::cout << std::endl;
}

