#pragma once

#include <eigen3/Eigen/Core>

namespace NuTo
{

//! @author Joerg F. Unger, ISM
//! @date September 2009
//! @brief abstract base class for the Transformations
class Transformation
{

public:
	virtual ~Transformation(){};
    //! @brief build the transformation using the given Points
	virtual void Build(const Eigen::MatrixXd& rCoordinates)=0;

    //! @brief transform the given points 
    virtual void TransformForward(Eigen::MatrixXd& rCoordinates)const=0;

    //! @brief transform the given points 
    virtual void TransformBackward(Eigen::MatrixXd& rCoordinates)const=0;

protected:

};

} // namespace nuto

