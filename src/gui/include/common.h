#ifndef __NUTOGUI_COMMON_H__
#define __NUTOGUI_COMMON_H__

// Provide a definition of nullptr if the compiler doesn't support it
#ifndef HAVE_NULLPTR
  #include "compat/nullptr"
#endif

// Macros for exporting/importing classes, symbols
#if defined(NUTOGUI_SHARED_LIBS)
  #define NUTOGUI_EXPORT	__attribute__((visibility("default")))
#endif

#ifndef NUTOGUI_EXPORT
  #define NUTOGUI_EXPORT
#endif
#ifndef NUTOGUI_IMPORT
  #define NUTOGUI_IMPORT
#endif

#endif // __NUTOGUI_COMMON_H__
