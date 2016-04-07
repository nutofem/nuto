# -*- coding: utf-8 -*-

import nuto

# Geometry/Mesh
area = 1.0
length = 1.0
num_elements = 10

# Material
conductivity = 1.0

# Boundaries
boundary_flux = 10.0
boundary_temperature = 20.0

# create one-dimensional structure
structure = nuto.Structure(1)

# create section
section = structure.SectionCreate("Truss")
structure.SectionSetArea(section, area)

# create material law
material = structure.ConstitutiveLawCreate("Heat_Conduction")
structure.ConstitutiveLawSetParameterDouble(material, "Thermal_Conductivity", conductivity)

# create nodes
node_coordinates = nuto.DoubleFullVector(1)
for node in range(0, num_elements + 1):
    node_coordinates.SetValue(0, 0, node * length/num_elements)
    structure.NodeCreate(node, node_coordinates)

# create interpolation type
interpolation_type = structure.InterpolationTypeCreate("Truss1D")
structure.InterpolationTypeAdd(interpolation_type, "coordinates", "equidistant1")
structure.InterpolationTypeAdd(interpolation_type, "temperature", "equidistant1")

# create elements
element_incidence = nuto.IntFullVector(2)
for element in range(0, num_elements):
    element_incidence.SetValue(0, 0, element)
    element_incidence.SetValue(1, 0, element + 1)
    structure.ElementCreate(interpolation_type, element_incidence)
    structure.ElementSetSection(element, section)
    structure.ElementSetConstitutiveLaw(element, material)

structure.ElementTotalConvertToInterpolationType()

# set boundary conditions and loads
direction = nuto.DoubleFullMatrix(1,1,(1,))
direction.SetValue(0, 0, 1.0);
structure.ConstraintLinearSetTemperatureNode(0, boundary_temperature)
structure.SetNumLoadCases(1)
structure.LoadCreateNodeHeatFlux(0, num_elements, direction, boundary_flux)

# start analysis
structure.SolveGlobalSystemStaticElastic(0)
intGradient = structure.BuildGlobalInternalGradient().J.Get("Temperature")
extGradientTmp = structure.BuildGlobalExternalLoadVector(0)
extGradient = extGradientTmp.J.Get("Temperature")
residual = nuto.DoubleFullVector(intGradient - extGradient)
print("Residual: {0}".format(residual.Norm()))

## visualize results
visualization_group = structure.GroupCreate("Elements");
structure.GroupAddElementsTotal(visualization_group)

structure.AddVisualizationComponent(visualization_group, "Temperature");
structure.ExportVtkDataFileElements("Temperature1D.vtk")
