#pragma once
#include<eigen3/Eigen/Dense>
#include "TypeDefs.h"

namespace NuTo
{

class CellInterpolationBase
{
public:
    CellInterpolationBase(){}

    //! @brief extracts all node values of this element
    //! @remark virtual to make it testable
    virtual NodeValues ExtractNodeValues() const = 0;

//    Eigen::VectorXd Interpolate(const NaturalCoords& rNaturalCoords) const
//    {
//        return mInterpolation.GetN(rNaturalCoords) * ExtractNodeValues();
//    }

//    NMatrix GetN(const NaturalCoords& rNaturalIPCoords) const
//    {
//        int dim = GetDofDimension();
//        Eigen::MatrixXd N(dim, dim * GetNumNodes());

//        auto shapeFunctions = GetShapeFunctions(rNaturalIPCoords);

//        for (int i = 0; i < GetNumNodes(); ++i)
//            N.block(0, i * dim, dim, dim) = Eigen::MatrixXd::Identity(dim, dim) * shapeFunctions[i];
//        return N;
//    }

    virtual int GetDofDimension() const = 0;

    virtual int GetNumNodes() const = 0;

    virtual Eigen::VectorXd GetShapeFunctions(Eigen::VectorXd ipCoords) const = 0;

    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(Eigen::VectorXd ipCoords) const = 0;
};
} /* NuTo */
