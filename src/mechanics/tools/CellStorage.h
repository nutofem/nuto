#pragma once
#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementCollection.h"
#include "base/Group.h"

namespace NuTo
{

class IntegrationTypeBase;

//! stores and creates integration cells
class CellStorage
{
public:
    //! creates and stores cells from elements, an integration type and a continuous id
    //! @param elements group of elements
    //! @param integrationType suitable integration type that matches the elements
    //! @param cellStartId start id of the continous cell numbering
    //! @return group of newly created cells
    Group<CellInterface> AddCells(Group<ElementCollectionFem> elements, const IntegrationTypeBase& integrationType,
                                  int cellStartId = 0);

private:
    boost::ptr_vector<CellInterface> mCells;
};
} /* NuTo */
