// $Id$ 
#ifndef SECTIONENUM_H_
#define SECTIONENUM_H_

namespace NuTo
{
namespace Section
{
enum eSectionType
{
    TRUSS,	          //!< one-dimensional truss behavior
    PLANE_STRAIN,     //!< two-dimensional behavior plane strain
    PLANE_STRESS,     //!< two-dimensional behavior plane stress
    VOLUME            //!< three-dimensional
};
}
}
#endif /* SECTIONENUM_H_ */
