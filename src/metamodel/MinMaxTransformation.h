#pragma once

#include "metamodel/Transformation.h"

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief abstract base class for the Transformations
class MinMaxTransformation : public Transformation
{

public:
    //! @brief constructor
    //! @param rCoordinate coordinate within the point coordinates (0<=entry<dim)
    //! @param rLb lower bound after the transformation
    //! @param rUb upper bound after the transformation
    MinMaxTransformation(unsigned int rCoordinate, double rLb, double rUb);

    //! @brief copy constructor
    //! @param rOther object which is copied
    MinMaxTransformation(const MinMaxTransformation& rOther);

    //! @brief destructor
    ~MinMaxTransformation()
    {
    }

    //! @brief build the transformation using the given Points
    //! @param rCoordinates input point coordinates
    virtual void Build(const Eigen::MatrixXd& rCoordinates) override;

    //! @brief transform the given points in forward direction x = f(x)
    //! @param rCoordinates input point coordinates
    virtual void TransformForward(Eigen::MatrixXd& rCoordinates) const override;

    //! @brief transform the given points in backward direction x = f^(-1)(x)
    //! @param rCoordinates input point coordinates
    virtual void TransformBackward(Eigen::MatrixXd& rCoordinates) const override;

protected:
    int mCoordinate; //!< coordinate within the point coordinates (0<=entry<dim
    double mMin; //!< min value of given coordinates
    double mMax; //!< max value of given coordinates
    double mUb; //!< upper bound after the transformation
    double mLb; //!< lower bound after the transformation

    // default constructor required by serialize
    MinMaxTransformation()
        : Transformation()
    {
    }
};
} // namespace nuto
