#pragma once
#include <eigen3/Eigen/Core>
#include "mechanics/interpolation/InterpolationEnum.h"

namespace NuTo
{

//! @brief N, the 'blown up' shape functions
typedef Eigen::MatrixXd NMatrix;

//! @brief BGradient in local coordinates. has to be transformed to BGradientGlobal using the jacobian
typedef Eigen::MatrixXd BGradientLocal;

//! @brief BStrain in local coordinates. has to be transformed to BStrainGlobal using the jacobian
typedef Eigen::MatrixXd BStrainLocal;


//! @brief Base class for the interpolation. The derived classes provide information about the actual interpolation.
//!        [TODO] Implement 'caching' of the results. Requests for the same local node coordinates
//!        must be calculated only once.
class Interpolation
{
public:
    Interpolation(eInterpolation rType, int rOrder, int rDofDimension)
        : mType(rType)
        , mOrder(rOrder)
        , mDofDimension(rDofDimension)
    {
    }
    virtual ~Interpolation()            = default; // virtual destructor needed, rule of 5 below...
    Interpolation(const Interpolation&) = default;
    Interpolation(Interpolation&&)      = default;
    Interpolation& operator=(const Interpolation&) = default;
    Interpolation& operator=(Interpolation&&) = default;

    int GetDofDimension() const
    {
        return mDofDimension;
    }

    NMatrix GetN(const Eigen::VectorXd& rLocalIPCoords) const;
    BGradientLocal GetBGradient(const Eigen::VectorXd& rLocalIPCoords) const;
    BStrainLocal GetBStrain(const Eigen::VectorXd& rLocalIPCoords) const;

    //! @brief calculates the shape functions
    //! @param rLocalIPCoords integration point coordinates in the local coordinate system
    //! @return vector of shape functions, dimension: [GetNumNodes() x 1]
    virtual Eigen::VectorXd GetShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const = 0;

    //! @brief calculates the derivative shape functions
    //! @param rLocalIPCoords integration point coordinates in the local coordinate system
    //! @return matrix of derivate shape functions, dimension: [GetNumNodes() x local dimension]
    //! @remark 'local dimension' above is the dimension of rLocalIPCoords
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(const Eigen::VectorXd& rLocalIPCoords) const = 0;

    //! @brief returns the local node coordinates
    //! @param rNodeId local node number
    //! @return node coordinates in the local coordinate system
    virtual Eigen::VectorXd GetLocalCoords(int rNodeId) const = 0;

    //! @brief returns the number of nodes
    //! @return number of nodes
    virtual int GetNumNodes() const = 0;


protected:
    eInterpolation mType;
    int mOrder;
    int mDofDimension;
};
} /* NuTo */
