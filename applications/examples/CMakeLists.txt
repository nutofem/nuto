unset(all_examples CACHE)

#add_subdirectory(c++)
#add_subdirectory(python)

foreach(test
        GradientDamageCalibration
    )

    add_executable(${test} EXCLUDE_FROM_ALL ${test}.cpp)
    target_link_libraries(${test} NuTo) 
    add_test(examples::cpp::${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
    append_to_examples(${test})
endforeach()

add_custom_target(examples DEPENDS ${all_examples})
