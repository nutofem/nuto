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

enum class eSubMatrix
{
    JJ,
    JK,
    KJ,
    KK
};

enum class eOutput
{
    DAMPING_FULL,
    DAMPING_MOISTURE_TRANSPORT,
    INTERNAL_GRADIENT,
    INTERNAL_GRADIENT_MOISTURE_TRANSPORT,
    MASS_FULL,
    STIFFNESS_FULL,
    STIFFNESS_DISPLACEMENTS_JJ,
    STIFFNESS_MOISTURE_TRANSPORT
};

}
}
#endif /* STRUCTUREBASEENUM_H_ */
