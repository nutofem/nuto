#pragma once

#include <sstream>

#include "nuto/mechanics/cell/CellInterface.h"
#include "nuto/mechanics/elements/ElementCollection.h"
#include "nuto/mechanics/dofs/DofContainer.h"
#include "nuto/mechanics/integrationtypes/IntegrationTypeBase.h"
#include "nuto/mechanics/cell/CellData.h"
#include "nuto/mechanics/cell/CellIpData.h"

namespace NuTo
{
class Cell : public CellInterface
{
public:
    Cell(const ElementCollection& elements, const IntegrationTypeBase& integrationType, const int id)
        : mElements(elements)
        , mIntegrationType(integrationType)
        , mId(id)
        , mShape(elements.GetShape())
    {
        if (elements.GetShape() != integrationType.GetShape())
        {
            std::stringstream message;
            message << "The shape of the element (" << elements.GetShape() << ") and integrationtype ("
                    << integrationType.GetShape() << ") don't match.";
            throw Exception(__PRETTY_FUNCTION__, message.str());
        }
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
            Jacobian jacobian(mElements.CoordinateElement().ExtractCoordinates(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(cellData, jacobian, ipCoords, iIP);
            f(cellipData);
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
        return NuTo::Interpolate(mElements.CoordinateElement(), naturalCoords);
    }

    Eigen::VectorXd Interpolate(Eigen::VectorXd naturalCoords, DofType dof) const override
    {
        return NuTo::Interpolate(mElements.DofElement(dof), naturalCoords);
    }

    std::vector<Eigen::VectorXd> Eval(EvalFunction f) const override
    {
        CellData cellData(mElements, Id());
        std::vector<Eigen::VectorXd> result;
        for (int iIP = 0; iIP < mIntegrationType.GetNumIntegrationPoints(); ++iIP)
        {
            auto ipCoords = mIntegrationType.GetLocalIntegrationPointCoordinates(iIP);
            Jacobian jacobian(mElements.CoordinateElement().ExtractCoordinates(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(cellData, jacobian, ipCoords, iIP);
            result.push_back(f(cellipData));
        }
        return result;
    }

    const Shape& GetShape() const override
    {
        return mShape;
    }

private:
    //! @brief integrates various operations with various return types
    //! @param f operation to perform
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
            Jacobian jacobian(mElements.CoordinateElement().ExtractCoordinates(),
                              mElements.CoordinateElement().GetDerivativeShapeFunctions(ipCoords));
            CellIpData cellipData(cellData, jacobian, ipCoords, iIP);
            result += f(cellipData) * jacobian.Det() * ipWeight;
        }
        return result;
    }

private:
    const ElementCollection& mElements;
    const IntegrationTypeBase& mIntegrationType;
    const int mId;
    const Shape& mShape;
};
} /* NuTo */
