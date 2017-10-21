#include "BoostUnitTest.h"
#include "mechanics/feti/ReverseMap.h"
#include "iostream"

BOOST_AUTO_TEST_CASE(Constructor)
{
    std::map<int, int> map;
    map.emplace(1, 11);
    map.emplace(2, 22);
    map.emplace(3, 22);

    NuTo::ReverseMap<int> reverseMap(map);

    std::cout << "Original map" << std::endl;
    for (const auto& pair : map)
    {
        std::cout << pair.first << "\t=>\t" << pair.second << std::endl;
    }

    std::cout << "Reverse map" << std::endl;
    for (const auto& pair : reverseMap)
    {
        std::cout << pair.first << "\t=>\t";
        for (const auto& i : pair.second)
            std::cout << i << ",\t";

        std::cout << std::endl;
    }

    std::cout << "reverseMap[22][0]: \t" << reverseMap[22][0] << std::endl;
    std::cout << "reverseMap[22][1]: \t" << reverseMap[22][1] << std::endl;

    BOOST_REQUIRE(map.size() == 3);
    BOOST_REQUIRE(reverseMap.size() == 2);
    BOOST_REQUIRE(reverseMap[11].size() == 1);
    BOOST_REQUIRE(reverseMap[22].size() == 2);
}


BOOST_AUTO_TEST_CASE(Add_map)
{

    std::map<int, int> map;
    map.emplace(1, 11);
    map.emplace(2, 22);
    map.emplace(3, 22);

    NuTo::ReverseMap<int> reverseMap(map);
    reverseMap.addMap(map);

    std::cout << "Reverse map after adding the original map a second time" << std::endl;
    for (const auto& pair : reverseMap)
    {
        std::cout << pair.first << "\t=>\t";
        for (const auto& i : pair.second)
            std::cout << i << ",\t";

        std::cout << std::endl;
    }

    BOOST_REQUIRE(map.size() == 3);
    BOOST_REQUIRE(reverseMap.size() == 2);
    BOOST_REQUIRE(reverseMap[11].size() == 2);
    BOOST_REQUIRE(reverseMap[22].size() == 4);
}
