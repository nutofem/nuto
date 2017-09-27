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

    DofVector<double> operator()(const VectorOperation& op) override
    {
        return Integrate(op, DofVector<double>());
    }

    DofMatrix<double> operator()(const MatrixOperation& op) override
    {
        return Integrate(op, DofMatrix<double>());
    }

    double operator()(const ScalarOperation& op) override
    {
        return Integrate(op, double{0});
    }

    void operator()(const VoidOperation& op) override
    {
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                              mCoordinateElement.GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(mElements, jacobian, ipCoords);
            op(mIntegrands[iIP], cellData, cellipData);
        }
    }

    DofVector<int> DofNumbering() override
    {
        return DofVector<int>();
    }

private:
    //! @brief integrates various operations with various return types
    //! @param op operation to perform
    //! @param result result value. It is not clear how to properly initialize an arbitrary TResult to zero. Thus, the
    //! user has to provide it with this argument.
    template <typename TOperation, typename TReturn>
    TReturn Integrate(TOperation&& op, TReturn result)
    {
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian jacobian(mCoordinateElement.ExtractNodeValues(),
                              mCoordinateElement.GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(mElements, jacobian, ipCoords);
            result += op(mIntegrands[iIP], cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return result;
    }

    const ElementInterface& mCoordinateElement;
    DofContainer<ElementInterface*> mElements;
    const IntegrationTypeBase& mIntegrationType;
    boost::ptr_vector<Integrand::Base> mIntegrands;
};
} /* NuTo */
