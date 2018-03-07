unset(all_examples CACHE)

#add_subdirectory(c++)
#add_subdirectory(python)
    
create_symlink(SingleInclusion.msh)
add_executable(ConcreteCompressionWithInclusion ConcreteCompressionWithInclusion.cpp)
target_link_libraries(ConcreteCompressionWithInclusion Base Math Mechanics Visualize Boost::filesystem)
append_to_examples(ConcreteCompressionWithInclusion)

add_custom_target(examples DEPENDS ${all_examples})
