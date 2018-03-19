#pragma once

#include <string>

namespace NuTo
{
namespace Visualize
{
class UnstructuredGrid;

class XMLWriter
{
public:
    static void Export(std::string filename, const UnstructuredGrid& unstructuredGrid, bool asBinary);
};
} /* Visualize */
} /* NuTo */
