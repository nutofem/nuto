#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

struct GmshNode;
struct GmshElement;

namespace NuTo
{


class MeshGmsh
{
public:
    explicit MeshGmsh(const std::string& fileName);

private:
    void ReadGmshFile(const std::string& fileName);
};
}
