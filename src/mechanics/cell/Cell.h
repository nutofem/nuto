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

    //! @brief builds the internal gradien
    DofVector Gradient()
    {
        DofVector gradient;
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                              mCoordinateElement.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));
            CellIPData cellipData(mElements, jacobian, ipCoords);
            gradient += mIntegrand->Gradient(cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return gradient;
    }

    //! @brief Extracts a vector (each IP) of vectors (several IPValues for the same integrion point) of IPValues
    std::vector<std::vector<IPValue>> IPValues()
    {
        std::vector<std::vector<IPValue>> ipValues;

        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                              mCoordinateElement.GetInterpolation().GetDerivativeShapeFunctions(ipCoords));
            CellIPData cellipData(mElements, jacobian, ipCoords);
            ipValues.push_back(mIntegrand->IPValues(cellData, cellipData));
        }
        return ipValues;
    }

private:

    const ElementSimple& mCoordinateElement;
    DofContainer<ElementSimple*> mElements;
    const IntegrationTypeBase& mIntegrationType;
    std::unique_ptr<Integrand> mIntegrand;
};
} /* NuTo */
