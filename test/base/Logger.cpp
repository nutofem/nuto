#include "BoostUnitTest.h"
#include "nuto/base/Logger.h"


std::string FromFile(std::string filename)
{
    std::string content;
    std::string line;
    std::ifstream inputfile(filename);
    while (getline(inputfile, line))
    {
        content += line;
        content += '\n';
    }
    return content;
}

BOOST_AUTO_TEST_CASE(Prefix)
{
    std::string filename = "prefixtest.log";
    NuTo::Log::Info.OpenFile(filename);
    NuTo::Log::Info.SetQuiet(true);

    NuTo::Log::Info << "Writing to file " << filename << ".\n";
    NuTo::Log::Info << "Writing again to file " << filename << '\n';
    NuTo::Log::Info << "Writing a \nmulti line log \nwith only one prefix.\n";

    std::string expected = R"([Info ] Writing to file prefixtest.log.
[Info ] Writing again to file prefixtest.log
[Info ] Writing a 
multi line log 
with only one prefix.)";

    BOOST_CHECK(FromFile(filename).compare(0, expected.size(), expected) == 0);
}
