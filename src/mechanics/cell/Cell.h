#pragma once

#include "mechanics/cell/CellInterface.h"
#include <boost/ptr_container/ptr_vector.hpp>
//#include "mechanics/elements/ElementSimple.h"
#include "mechanics/nodes/DofContainer.h"
#include "mechanics/cell/Integrand.h"
#include "mechanics/cell/PDE_Element.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/cell/PDE_Element.h"

namespace NuTo
{
template <class TIntegrand>
class Cell : public CellInterface
{
public:
    Cell(const PDE_Element& element, const IntegrationTypeBase& integrationType, const TIntegrand& integrand)
        : mElement(element)
        , mIntegrationType(integrationType)
        , mIntegrandVec()
    {
        for (int i = 0; i < integrationType.GetNumIntegrationPoints(); i++)
            mIntegrandVec.push_back(integrand.Clone());
    }

    //! @brief builds the internal gradien
    DofVector<double> Gradient() override
    {
        DofVector<double> gradient;
        CellData cellData(mElement);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            CellIPData cellipData(mElement, ipCoords);
            gradient += mIntegrandVec[iIP].Gradient(cellData, cellipData) * cellipData.Jacobian().Det() * ipWeight;
        }
        return gradient;
    }

    //! @brief builds the hessian0 matrix
    DofMatrix<double> Hessian0() override
    {
        DofMatrix<double> hessian0;
        CellData cellData(mElement);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            CellIPData cellipData(mElement, ipCoords);
            hessian0 += mIntegrandVec[iIP].Hessian0(cellData, cellipData) * cellipData.Jacobian().Det() * ipWeight;
        }
        return hessian0;
    }


    DofVector<int> DofNumbering() override
    {
        return DofVector<int>();
    }

    //! @brief Extracts a vector (each IP) of vectors (several IPValues for the same integrion point) of IPValues
    std::vector<std::vector<IPValue>> IPValues() override
    {
        std::vector<std::vector<IPValue>> ipValues;

        CellData cellData(mElement);
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
             CellIPData cellipData(mElement, ipCoords);
            //ipValues.push_back(mIntegrand[iIP].IPValues(cellData, cellipData));
        }
        return ipValues;
    }

private:
    const PDE_Element& mElement;
    const IntegrationTypeBase& mIntegrationType;
    boost::ptr_vector<TIntegrand> mIntegrandVec;
};
} /* NuTo */
