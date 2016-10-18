#pragma once

namespace NuTo
{
enum class eSectionType
{
    TRUSS,	            //!< one-dimensional truss behavior
    PLANE_STRAIN,       //!< two-dimensional behavior plane strain
    PLANE_STRESS,       //!< two-dimensional behavior plane stress
    VOLUME,             //!< three-dimensional
    SPRING,             //!< one-dimensional spring
    FIBRE_MATRIX_BOND,  //!< section for the interface bewteen fibre and matrix material
};
}// namespace NuTo
