#include <iostream>
#include "nuto/math/MortonOrder.h"

int main()
{
    uint32_t x = 2;
    uint32_t y = 1;
    uint32_t z = 1;

    uint32_t key = NuTo::MortonOrder::EncodeMorton3D(x, y, z);
    std::cout << " result " << key << "\n";
    return 0;
}
