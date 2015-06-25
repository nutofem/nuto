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

enum class eOutput
{
    STIFFNESS_DISPLACEMENTS_JJ
};

}
}
#endif /* STRUCTUREBASEENUM_H_ */
