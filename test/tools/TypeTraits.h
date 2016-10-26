//
// Created by Thomas Titscher on 10/26/16.
//
#pragma once
#include <type_traits>
namespace NuTo
{
namespace Test
{
template <typename T>
struct Copy
{
    static_assert(std::is_copy_constructible<T>::value, "The class is not copy constructable.");
    static_assert(std::is_copy_assignable<T>::value, "The class is not copy assignable.");
};
template <typename T>
struct Move
{
    static_assert(std::is_move_constructible<T>::value, "The class is not move constructable.");
    static_assert(std::is_move_assignable<T>::value, "The class is not move assignable.");
};
}
}
