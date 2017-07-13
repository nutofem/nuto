#pragma once

namespace NuTo
{
enum class eDirection
{
    X = 0,
    Y = 1,
    Z = 2
};

inline int ToComponentIndex(eDirection direction)
{
    return static_cast<int>(direction);
}

} // namespace NuTo
