#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/TypeDefs.h"

namespace NuTo
{

class ElementInterface
{
public:
    virtual ~ElementInterface() noexcept = default;

    //! @brief extracts all node values of this element
    //!
    //! They are ordered in blocks belonging to a node, e.g:
    //! coordinate Nodes in 3D:
    //!   NodeValues: (x1,y1,z1, x2,y2,z2, ...)
    virtual Eigen::VectorXd ExtractNodeValues(int instance = 0) const = 0;

    virtual int GetDofDimension() const = 0;

    //! @brief extract the dof numbers from its nodes.
    //! @remark They have to be in the same order as defined in ExtractNodeValues()
    virtual Eigen::VectorXi GetDofNumbering() const = 0;
    virtual int GetNumNodes() const = 0;
    virtual Eigen::MatrixXd GetNMatrix(NaturalCoords ipCoords) const = 0;
    virtual Eigen::VectorXd GetShapeFunctions(NaturalCoords ipCoords) const = 0;
    virtual Eigen::MatrixXd GetDerivativeShapeFunctions(NaturalCoords ipCoords) const = 0;
};

inline Eigen::VectorXd Interpolate(const ElementInterface& element, NaturalCoords ipCoords)
{
    return element.GetNMatrix(ipCoords) * element.ExtractNodeValues();
}


} /* NuTo */
