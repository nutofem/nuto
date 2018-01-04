#pragma once

#include "mechanics/cell/CellInterface.h"
#include "mechanics/elements/ElementCollection.h"
#include "mechanics/dofs/DofContainer.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/cell/CellData.h"
#include "mechanics/cell/CellIpData.h"

namespace NuTo
{
class Cell : public CellInterface
{
public:
    Cell(const ElementCollection& elements, const IntegrationTypeBase& integrationType, const int id)
        : mElements(elements)
        , mIntegrationType(integrationType)
        , mId(id)
    {
    }

    DofVector<double> Integrate(VectorFunction f) override
    {
        return IntegrateGeneric(f, DofVector<double>());
    }

    DofMatrix<double> Integrate(MatrixFunction f) override
    {
        return IntegrateGeneric(f, DofMatrix<double>());
    }

    double Integrate(ScalarFunction f) override
    {
        return IntegrateGeneric(f, double{0});
    }

    void Apply(VoidFunction f) override
    {
        CellData cellData(mElements, Id());
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mElements.CoordinateElement().ExtractNodeValues(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords),
                              mElements.CoordinateElement().GetDofDimension());
            CellIpData cellipData(mElements, jacobian, ipCoords, iIP);
            f(cellData, cellipData);
        }
    }

    Eigen::VectorXi DofNumbering(DofType dof) override
    {
        return mElements.DofElement(dof).GetDofNumbering();
    }

    int Id() const
    {
        return mId;
    }

    Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords) const override
    {
        const auto& coordinateElement = mElements.CoordinateElement();
        return coordinateElement.GetNMatrix(naturalCoords) * coordinateElement.ExtractNodeValues();
    }

    Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords, DofType dof) const override
    {
        const auto& element = mElements.DofElement(dof);
        return element.GetNMatrix(naturalCoords) * element.ExtractNodeValues();
    }

    std::vector<Eigen::VectorXd>
    Eval(std::function<Eigen::VectorXd(const CellData& cellData, const CellIpData& cellIpData)> f) const override
    {
        CellData cellData(mElements, Id());
        std::vector<Eigen::VectorXd> result;
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mElements.CoordinateElement().ExtractNodeValues(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords),
                              mElements.CoordinateElement().GetDofDimension());
            CellIpData cellipData(mElements, jacobian, ipCoords, iIP);
            result.push_back(f(cellData, cellipData));
        }
        return result;
    }

private:
    //! @brief integrates various operations with various return types
    //! @param op operation to perform
    //! @param result result value. It is not clear how to properly initialize an arbitrary TResult to zero. Thus, the
    //! user has to provide it with this argument.
    template <typename TOperation, typename TReturn>
    TReturn IntegrateGeneric(TOperation&& f, TReturn result)
    {
        CellData cellData(mElements, Id());
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            auto ipWeight = mIntegrationType.GetIntegrationPointWeight(iIP);
            Jacobian jacobian(mElements.CoordinateElement().ExtractNodeValues(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords),
                              mElements.CoordinateElement().GetDofDimension());
            CellIpData cellipData(mElements, jacobian, ipCoords, iIP);
            result += f(cellData, cellipData) * jacobian.Det() * ipWeight;
        }
        return result;
    }

private:
    const ElementCollection& mElements;
    const IntegrationTypeBase& mIntegrationType;
    const int mId;
};
} /* NuTo */
