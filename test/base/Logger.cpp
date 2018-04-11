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
    NuTo::Logger l("[Prefix]");

    std::string filename = "prefixtest.log";
    l.OpenFile(filename);
    l.SetQuiet(true);

    l << "Writing to file " << filename << ".\n";
    l << "Writing again to file " << filename << '\n';
    l << "Writing a \nmulti line log \nwith only one prefix.\n";

    std::string expected = R"([Prefix] Writing to file prefixtest.log.
[Prefix] Writing again to file prefixtest.log
[Prefix] Writing a 
multi line log 
with only one prefix.)";

    BOOST_CHECK(FromFile(filename).compare(0, expected.size(), expected) == 0);
}
