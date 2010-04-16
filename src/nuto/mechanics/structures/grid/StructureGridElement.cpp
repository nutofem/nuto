#include "nuto/mechanics/structures/grid/StructureGrid.h"
#include <sstream>

//! @brief returns the number of nodes
//! @return number of nodes
int NuTo::StructureGrid::GetNumElements() const
{
    return mElementVec.size();
}

//! @brief a reference to a node
//! @param identifier
//! @return reference to a node
NuTo::ElementBase* NuTo::StructureGrid::ElementGetElementPtr(int rElementNumber)
{
    if (rElementNumber<0 || rElementNumber>=GetNumElements())
        throw MechanicsException("[NuTo::StructureGrid::ElementGetElementPtr] Conversion from string to int did not yield valid element number.");
    return &mElementVec[rElementNumber];
}

//! @brief gives the identifier of a node
//! @param reference to a node
//! @return identifier
int NuTo::StructureGrid::ElementGetId(ElementBase* rElement)const
{
    int elementNumber(0);
    boost::ptr_vector<ElementBase>::const_iterator it;
    for (it = mElementVec.begin(); it!= mElementVec.end(); it++,elementNumber++)
    {
        if (&(*it)==rElement)
        {
            break;
        }
    }
    if (it==mElementVec.end())
        throw MechanicsException("[NuTo::StructureGrid::GetElementId] Element does not exist.");
    return elementNumber;
}
