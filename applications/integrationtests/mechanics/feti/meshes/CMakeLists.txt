################################################################################
##  MESH FILES
################################################################################


foreach(rank RANGE 0 3)
    foreach(mesh FetiMultiplicityScaling.mesh)
        create_symlink(${mesh}${rank})
    endforeach()
endforeach()

foreach(rank RANGE 0 1)
    foreach(mesh Feti.mesh)
        create_symlink(${mesh}${rank})
    endforeach()
endforeach()

foreach(rank RANGE 0 15)
    foreach(mesh FetiRixenFarhatHeterogeneousExample.mesh)
        create_symlink(${mesh}${rank})
    endforeach()
endforeach()

create_symlink(Feti_compare.msh)
