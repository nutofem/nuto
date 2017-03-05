#pragma once

#include <vector>
#include "mechanics/elements/Element.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/dofSubMatrixStorage/BlockFullVector.h"
#include "mechanics/cell/Integrand.h"

namespace NuTo
{
class Cell
{
public:
    Cell(const Element& rCoordinateElement, DofContainer<Element*> rElements,
          const Integrand& rIntegrand)
        : mCoordinateElement(rCoordinateElement)
        , mElements(rElements)
        , mIntegrand(rIntegrand.Clone())
    {
    }

    BlockFullVector<double> Gradient()
    {
    }


private:
    const Element& mCoordinateElement;
    DofContainer<Element*> mElements;
    std::unique_ptr<Integrand> mIntegrand;
};
} /* NuTo */
