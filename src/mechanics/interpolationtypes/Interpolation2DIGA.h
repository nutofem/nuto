#pragma once

#include "mechanics/interpolationtypes/InterpolationBaseIGA.h"

namespace NuTo
{

class Interpolation2DIGA: public InterpolationBaseIGA
{
#ifdef ENABLE_SERIALIZATION
    friend class boost::serialization::access;
    //! @brief default constructor for serialization
protected:
    Interpolation2DIGA(){}
#endif  // ENABLE_SERIALIZATION

public:

    Interpolation2DIGA(NuTo::Node::eDof rDofType,
                       NuTo::Interpolation::eTypeOrder rTypeOrder,
                       int rDimension,
                       const Eigen::Vector2i &rDegree,
                       const Eigen::VectorXd &rKnotsX,
                       const Eigen::VectorXd &rKnotsY,
                       const Eigen::MatrixXd &rWeights);

    int GetSplineDegree(int dir) const override
    {
        assert(dir == 0 || dir == 1);
        return mDegree(dir);
    }

    eIntegrationType GetStandardIntegrationType() const override;

    std::vector<Eigen::VectorXd> GetSurfaceEdgesCoordinates(int rSurface) const override;

    int GetLocalDimension() const override
    {
        return 2;
    }

    void UpdateIntegrationType(const IntegrationTypeBase& rIntegrationType) override;

    //********************************************
    //       SHAPE FUNCTIONS AND DERIVATIVES
    //********************************************

    // --- shape functions --- //

    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates) const override;

    //! @brief returns specific shape functions at a parameter, whicg fits the knot vector
    //! @param rCoordinates ... parameter
    //! @param rKnotIDs ... knot ids specifying the knot interval the rCoordinates are lying in (no need to search)
    //! @return ... specific shape functions
    Eigen::VectorXd CalculateShapeFunctions(const Eigen::VectorXd& rCoordinates, const Eigen::Vector2i &rKnotIDs) const;

    Eigen::VectorXd ShapeFunctionsIGA(int rIP, const Eigen::VectorXi &rKnotIDs) const override;

    // --- derivatives shape functions --- //

    Eigen::MatrixXd CalculateDerivativeShapeFunctionsNatural(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd DerivativeShapeFunctionsNaturalIGA(const Eigen::VectorXd &rCoordinates, const Eigen::VectorXi &rKnotIDs) const override;

    Eigen::MatrixXd DerivativeShapeFunctionsNaturalIGA(int rIP, const Eigen::VectorXi &rKnotIDs) const override;

    // --- N-matrix --- //

    Eigen::MatrixXd CalculateMatrixN(const Eigen::VectorXd& rCoordinates) const override;

    Eigen::MatrixXd MatrixNIGA(const Eigen::VectorXd& rCoordinates, const Eigen::VectorXi &rKnotIDs) const override;

    Eigen::MatrixXd MatrixNIGA(int rIP, const Eigen::VectorXi &rKnotIDs) const override;

    Eigen::MatrixXd MatrixNDerivativeIGA(const Eigen::VectorXd& rParameters, const Eigen::VectorXi& rKnotIDs, int rDerivative, int rDirection) const override;

    Eigen::MatrixXd ConstructMatrixN(Eigen::VectorXd rShapeFunctions) const;

    //********************************************
    //       SURFACE PARAMETRIZATION
    //********************************************

    Eigen::VectorXd CalculateNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface, const Eigen::MatrixXd &rKnots) const override;

    Eigen::MatrixXd CalculateDerivativeNaturalSurfaceCoordinates(const Eigen::VectorXd& rNaturalSurfaceCoordinates, int rSurface) const override;

    Eigen::VectorXi GetSurfaceNodeIndices(int rSurface) const override;

    int GetSurfaceDegree(int rSurface) const override;

    //! @brief returns the number of surfaces
    inline int GetNumSurfaces() const override
    {
        return 4;
    }


#ifdef ENABLE_SERIALIZATION
    //! @brief serializes the class
    //! @param ar         archive
    //! @param version    version
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
#endif  // ENABLE_SERIALIZATION

protected:

    //********************************************
    //               MEMBERS
    //********************************************

    //! @brief polynomial degree
    Eigen::Vector2i mDegree;

    //! @brief Knot vector x
    Eigen::VectorXd mKnotsX;

    //! @brief Knot vector y
    Eigen::VectorXd mKnotsY;

    //! @brief weights
    Eigen::MatrixXd mWeights;

    //! @brief integration point coordinates
    Eigen::Matrix<double, Eigen::Dynamic, 2> mIPCoordinates;

    //! @brief return the number node depending the shape and the order
    int CalculateNumNodes() const override;

};

} /* namespace NuTo */

#ifdef ENABLE_SERIALIZATION
BOOST_CLASS_EXPORT_KEY(NuTo::Interpolation2DIGA)
#endif
