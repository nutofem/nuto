#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/integrands/Base.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"

namespace NuTo
{
class Cell : public CellInterface
{
public:
    Cell(const ElementCollection& elements, const IntegrationTypeBase& integrationType,
         const Integrands::Base& integrand)
        : mElements(elements)
        , mIntegrationType(integrationType)
        , mIntegrands()
    {
        for (int i = 0; i < integrationType.GetNumIntegrationPoints(); i++)
            mIntegrands.push_back(integrand.Clone().release());
    }

    DofVector<double> Integrate(const VectorOperation& op) override
    {
        return Integrate(op, DofVector<double>());
    }

    DofMatrix<double> Integrate(const MatrixOperation& op) override
    {
        return Integrate(op, DofMatrix<double>());
    }

    double Integrate(const ScalarOperation& op) override
    {
        return Integrate(op, double{0});
    }

    void Apply(const VoidOperation& op) override
    {
        CellData cellData(mElements);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mElements.CoordinateElement().ExtractNodeValues(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords),
                              mElements.CoordinateElement().GetDofDimension());
            CellIpData cellipData(mElements, jacobian, ipCoords);
            op(mIntegrands[iIP], cellData, cellipData);
        }
    }

    Eigen::VectorXi DofNumbering(DofType dof) override
    {
        return mElements.DofElement(dof).GetDofNumbering();
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
            Jacobian jacobian(mElements.CoordinateElement().ExtractNodeValues(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords),
                              mElements.CoordinateElement().GetDofDimension());
            CellIpData cellipData(mElements, jacobian, ipCoords);
            result += op(mIntegrands[iIP], cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return result;
    }

private:
    const ElementCollection& mElements;
    const IntegrationTypeBase& mIntegrationType;
    boost::ptr_vector<Integrands::Base> mIntegrands;
};
} /* NuTo */
