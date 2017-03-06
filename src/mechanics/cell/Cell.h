#pragma once

#include <vector>
#include "mechanics/elements/ElementSimple.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/nodes/DofVector.h"
#include "mechanics/cell/Integrand.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"

namespace NuTo
{
class Cell
{
public:
    Cell(const ElementSimple& rCoordinateElement, DofContainer<ElementSimple*> rElements,
         const IntegrationTypeBase& rIntegrationType, const Integrand& rIntegrand)
        : mCoordinateElement(rCoordinateElement)
        , mElements(rElements)
        , mIntegrationType(rIntegrationType)
        , mIntegrand(rIntegrand.Clone())
    {
    }

    DofVector<Eigen::VectorXd> Gradient()
    {
        auto ipCoordinate = mIntegrationType.GetLocalIntegrationPointCoordinates(0);
        DofVector<Eigen::VectorXd> gradient;
        return gradient;
    }


private:
    const ElementSimple& mCoordinateElement;
    DofContainer<ElementSimple*> mElements;
    const IntegrationTypeBase& mIntegrationType;
    std::unique_ptr<Integrand> mIntegrand;
};
} /* NuTo */
