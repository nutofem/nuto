#include "mechanics/tools/CellStorage.h"
#include "mechanics/cell/Cell.h"

using namespace NuTo;

Group<CellInterface> CellStorage::AddCells(Group<ElementCollectionFem> elements,
                                           const IntegrationTypeBase& integrationType, int cellStartId)
{
    Group<CellInterface> cellGroup;
    for (auto& element : elements)
    {
        mCells.push_back(new Cell(element, integrationType, cellStartId++));
        cellGroup.Add(mCells.back());
    }
    return cellGroup;
}
