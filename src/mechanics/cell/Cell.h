#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementInterface.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/integrands/Base.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"

namespace NuTo
{
class Cell : public CellInterface
{
public:
    Cell(const ElementInterface& coordinateElement, DofContainer<ElementInterface*> elements,
         const IntegrationTypeBase& integrationType, const Integrand::Base& integrand)
        : mCoordinateElement(coordinateElement)
        , mElements(elements)
        , mIntegrationType(integrationType)
        , mIntegrands()
    {
        for (int i = 0; i < integrationType.GetNumIntegrationPoints(); i++)
            mIntegrands.push_back(integrand.Clone().release());
    }

    DofVector<double> Integrate(const VectorOperation& op) override
    {
        DofVector<double> gradient;
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                                    mCoordinateElement.GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(mElements, jacobian, ipCoords);
            gradient += op(mIntegrands[iIP], cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return gradient;
    }

    DofMatrix<double> Integrate(const MatrixOperation& op) override
    {
        DofMatrix<double> hessian0;
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                                    mCoordinateElement.GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(mElements, jacobian, ipCoords);
            hessian0 += op(mIntegrands[iIP], cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return hessian0;
    }

    DofVector<int> DofNumbering() override
    {
        return DofVector<int>();
    }

private:
    const ElementInterface& mCoordinateElement;
    DofContainer<ElementInterface*> mElements;
    const IntegrationTypeBase& mIntegrationType;
    boost::ptr_vector<Integrand::Base> mIntegrands;
};
} /* NuTo */
