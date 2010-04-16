#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include "nuto/math/FullMatrix.h"
#include <sstream>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumNodes() const
{
    return mNodeVec.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::NodeBase* NuTo::StructureGrid::NodeGetNodePtr(int rNodeNumber)
{
    if (rNodeNumber<0 || rNodeNumber>=GetNumNodes())
        throw MechanicsException("[NuTo::StructureGrid::NodeGetNodePtr] Conversion from string to int did not yield valid node number.");
    return &mNodeVec[rNodeNumber];
}

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::NodeGetId(NodeBase* rNode)const
{
    int nodeNumber(0);
    boost::ptr_vector<NodeBase>::const_iterator it;
    for (it = mNodeVec.begin(); it!= mNodeVec.end(); it++,nodeNumber++)
    {
        if (&(*it)==rNode)
        {
            break;
        }
    }
    if (it== mNodeVec.end())
        throw MechanicsException("[NuTo::StructureGrid::GetNodeId] Node does not exist.");
    return nodeNumber;
}

// extract dof values (e.g. displacements, temperatures to the nodes)
void NuTo::StructureGrid::NodeExtractDofValues(FullMatrix<double>& rActiveDofValues, FullMatrix<double>& rDependentDofValues) const
{
    if (this->mNodeNumberingRequired)
    {
        throw MechanicsException("[NuTo::GridStructure::NodeExtractDofValues] a valid dof numbering was not found (build dof numbering using NodeBuildGlobalDofs).");
    }
    rActiveDofValues.Resize(this->mNumActiveDofs,1);
    rDependentDofValues.Resize(this->mNumDofs - this->mNumActiveDofs,1);

    // extract dof values from nodes
    for (boost::ptr_vector<NodeBase>::const_iterator it = this->mNodeVec.begin(); it!= this->mNodeVec.end(); it++)
    {
        it->GetGlobalDofValues(rActiveDofValues, rDependentDofValues);
    }
}
