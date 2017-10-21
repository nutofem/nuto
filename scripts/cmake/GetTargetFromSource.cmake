# in three steps go from "/home/user/nuto/src/mechanics/MyClass.cpp"
# to "mechanics.MyClass", and export this to Target
function(get_target_from_source Source Target)
    string(REPLACE "${CMAKE_SOURCE_DIR}/src/" "" libName ${Source})
    string(REPLACE ".cpp" "" libName ${libName})
    string(REPLACE "/" "." libName ${libName})
    set(${Target} ${libName} PARENT_SCOPE)
endfunction()
