// $Id$ 
#ifndef NODEENUM_H_
#define NODEENUM_H_

namespace NuTo
{
namespace Node
{
enum eNodeType
{
    NodeCoordinates,
    NodeCoordinatesDof,
    NodeCoordinatesMultiscale2D,
    NodeCoordinatesDofNonlocalData,
    NodeCoordinatesDofRadius,
    NodeGrid1D,
    NodeGrid2D,
    NodeGrid3D,
    NodeGridDisplacements1D,
    NodeGridDisplacements2D,
    NodeGridDisplacements3D
};

enum eAttributes
{
    COORDINATES,
    ROTATIONS,
    TEMPERATURES,
    DISPLACEMENTS,
    FINESCALEDISPLACEMENTS,
    NONLOCALDATA,
    NONLOCALEQPLASTICSTRAIN,
    NONLOCALTOTALSTRAIN,
    NONLOCALEQSTRAIN,
    WATERPHASEFRACTION,
    RELATIVEHUMIDITY
};

}
}
#endif /* NODEENUM_H_ */
