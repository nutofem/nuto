#pragma once
namespace NuTo
{
struct Voigt
{
    static constexpr int Dim(int d)
    {
        return d == 1 ? 1 : (d == 2 ? 3 : 6);
    }
};
} /* NuTo */
