import nuto
import numpy as np

structure1D = nuto.Structure(1)
structure2D = nuto.Structure(2)
structure3D = nuto.Structure(3)

# create section
area = 1.0
thickness = 1.0
section1D = structure1D.SectionCreate("Truss")
section2D = structure2D.SectionCreate("Plane_Strain")
section3D = structure3D.SectionCreate("Volume")
structure1D.SectionSetArea(section1D, area)
structure2D.SectionSetThickness(section2D, thickness)

# create material law
conductivity = 1.0
material1D = structure1D.ConstitutiveLawCreate("Heat_Conduction")
material2D = structure2D.ConstitutiveLawCreate("Heat_Conduction")
material3D = structure3D.ConstitutiveLawCreate("Heat_Conduction")
structure1D.ConstitutiveLawSetParameterDouble(material1D, "Thermal_Conductivity", conductivity)
structure2D.ConstitutiveLawSetParameterDouble(material2D, "Thermal_Conductivity", conductivity)
structure3D.ConstitutiveLawSetParameterDouble(material3D, "Thermal_Conductivity", conductivity)

# create nodes
node_id = 0
for x in range(2):
    structure1D.NodeCreate(node_id, np.array([float(x)]))
    node_id += 1

node_id = 0
for x in range(2):
    for y in range(2):
        structure2D.NodeCreate(node_id, np.array([float(x), float(y)]))
        node_id += 1

node_id = 0
for x in range(2):
    for y in range(2):
        for z in range(2):
            structure3D.NodeCreate(node_id, np.array([float(x), float(y), float(z)]))
            node_id += 1

# create interpolation type
truss_iptype = structure1D.InterpolationTypeCreate("Truss1D")
structure1D.InterpolationTypeAdd(truss_iptype, "coordinates", "equidistant1")
structure1D.InterpolationTypeAdd(truss_iptype, "temperature", "equidistant1")

quad_iptype = structure2D.InterpolationTypeCreate("Quad2D")
structure2D.InterpolationTypeAdd(quad_iptype, "coordinates", "equidistant1")
structure2D.InterpolationTypeAdd(quad_iptype, "temperature", "equidistant1")

brick_iptype = structure3D.InterpolationTypeCreate("Brick3D")
structure3D.InterpolationTypeAdd(brick_iptype, "coordinates", "equidistant1")
structure3D.InterpolationTypeAdd(brick_iptype, "temperature", "equidistant1")

# create one truss element
element_incidence = nuto.IntVector(2)
element_incidence[0] = 0
element_incidence[1] = 1
structure1D.ElementCreate(truss_iptype, element_incidence)
structure1D.ElementSetSection(0, section1D)
structure1D.ElementSetConstitutiveLaw(0, material1D)

# create one quad element
element_incidence = nuto.IntVector(4)
element_incidence[0] = 0
element_incidence[1] = 1
element_incidence[2] = 3
element_incidence[3] = 2
structure2D.ElementCreate(quad_iptype, element_incidence)
structure2D.ElementSetSection(0, section2D)
structure2D.ElementSetConstitutiveLaw(0, material2D)

# create one brick element
element_incidence = nuto.IntVector(8)
element_incidence[0] = 0
element_incidence[1] = 1
element_incidence[2] = 3
element_incidence[3] = 2
element_incidence[4] = 4
element_incidence[5] = 5
element_incidence[6] = 7
element_incidence[7] = 6
structure3D.ElementCreate(brick_iptype, element_incidence)
structure3D.ElementSetSection(0, section3D)
structure3D.ElementSetConstitutiveLaw(0, material3D)

structure1D.ElementTotalConvertToInterpolationType()
structure2D.ElementTotalConvertToInterpolationType()
structure3D.ElementTotalConvertToInterpolationType()

# check stiffness matrices
delta = 1e-6;
rel_tolerance = 1e-4;
structure1D.CheckHessian0(delta, rel_tolerance, True)
structure2D.CheckHessian0(delta, rel_tolerance, True)
structure3D.CheckHessian0(delta, rel_tolerance, True)




