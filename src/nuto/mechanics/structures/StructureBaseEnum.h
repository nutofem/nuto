// $Id$
#ifndef STRUCTUREBASEENUM_H_
#define STRUCTUREBASEENUM_H_

namespace NuTo
{
namespace StructureEnum
{

enum class eMatrixType
{
    STIFFNESS=0,
    DAMPING,
    MASS,
    LUMPED_MASS
};

enum class eSubVector
{
    J,
    K
};

enum class eSubMatrix
{
    JJ,
    JK,
    KJ,
    KK
};

enum class eOutput
{
    DAMPING,
    INTERNAL_GRADIENT,
    MASS,
    RESIDUAL_NORM_FACTOR,
    STIFFNESS,
    STIFFNESS_DISPLACEMENTS_JJ,
};

}
}
#endif /* STRUCTUREBASEENUM_H_ */
