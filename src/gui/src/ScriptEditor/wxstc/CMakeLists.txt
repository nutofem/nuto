SET(stc_wxSTC_SourceDir "${CMAKE_CURRENT_SOURCE_DIR}/src/stc")
SET(stc_Scintilla_SourceDir "${stc_wxSTC_SourceDir}/scintilla")

set_directory_properties(PROPERTIES
  INCLUDE_DIRECTORIES "${CMAKE_CURRENT_SOURCE_DIR}/include;${wxWidgets_INCLUDE_DIRS};${stc_Scintilla_SourceDir}/include;${stc_Scintilla_SourceDir}/src")

file(GLOB stc_wxSTC_sources "${stc_wxSTC_SourceDir}/*.h" "${stc_wxSTC_SourceDir}/*.cpp")
file(GLOB stc_Scintilla_sources "${stc_Scintilla_SourceDir}/src/*.h" "${stc_Scintilla_SourceDir}/src/*.cxx")

add_library(wxstc ${stc_wxSTC_sources} ${stc_Scintilla_sources})
set_target_properties(wxstc PROPERTIES COMPILE_DEFINITIONS "SCI_LEXER;LINK_LEXERS;__WX__")
