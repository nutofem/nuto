#pragma once
#include <eigen3/Eigen/Core>
#include "nuto/mechanics/interpolation/TypeDefs.h"
#include "nuto/mechanics/elements/ElementInterface.h"

namespace NuTo
{

class CoordinateElementInterface : public ElementInterface
{
public:
    virtual ~CoordinateElementInterface() noexcept = default;

    //! @brief extracts all node values of this element
    //!
    //! They are ordered in blocks belonging to a node, e.g:
    //! coordinate Nodes in 3D:
    //!   NodeValues: (x1,y1,z1, x2,y2,z2, ...)
    virtual Eigen::VectorXd ExtractCoordinates() const = 0;
};

inline Eigen::VectorXd Interpolate(const CoordinateElementInterface& element, NaturalCoords ipCoords)
{
    return element.GetNMatrix(ipCoords) * element.ExtractCoordinates();
}


} /* NuTo */
