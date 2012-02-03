// $Id$ 
#ifndef NODEENUM_H_
#define NODEENUM_H_

namespace NuTo
{
namespace Node
{
enum eNodeType
{
    NodeAccelerations1D,
    NodeAccelerations2D,
    NodeAccelerations3D,
    NodeCoordinates1D,
    NodeCoordinates2D,
    NodeCoordinates3D,
    NodeCoordinatesDisplacements1D,
    NodeCoordinatesDisplacements2D,
    NodeCoordinatesDisplacements3D,
    NodeCoordinatesDisplacementsMultiscale1D,
    NodeCoordinatesDisplacementsMultiscale2D,
    NodeCoordinatesDisplacementsMultiscale3D,
    NodeCoordinatesDisplacementsNonlocalData2D,
    NodeCoordinatesDisplacementsNonlocalData3D,
    NodeCoordinatesDisplacementsRadius2D,
    NodeCoordinatesDisplacementsRadius3D,
    NodeCoordinatesDisplacementsRotations2D,
    NodeCoordinatesDisplacementsRotations3D,
    NodeCoordinatesDisplacementsRotationsRadius2D,
    NodeCoordinatesDisplacementsRotationsRadius3D,
    NodeTemperature,
    NodeCoordinatesDisplacementsVelocitiesAccelerations1D,
    NodeCoordinatesDisplacementsVelocitiesAccelerations2D,
    NodeCoordinatesDisplacementsVelocitiesAccelerations3D,
    NodeCoordinatesTemperature1D,
    NodeCoordinatesTemperature2D,
    NodeCoordinatesTemperature3D,
    NodeDisplacements1D,
    NodeDisplacements2D,
    NodeDisplacements3D,
    NodeDisplacementsMultiscale2D,
    NodeGrid1D,
    NodeGrid2D,
    NodeGrid3D,
    NodeGridDisplacements1D,
    NodeGridDisplacements2D,
    NodeGridDisplacements3D,
    NodeRadius,
    NodeRotations2D,
    NodeRotations3D,
    NodeVelocities1D,
    NodeVelocities2D,
    NodeVelocities3D
};

enum eAttributes
{
    ACCELERATIONS=0,
    COORDINATES,
    DISPLACEMENTS,
    NONLOCALDATA,
    RADIUS,
    ROTATIONS,
    TEMPERATURES,
    VELOCITIES
};
}
}
#endif /* NODEENUM_H_ */
