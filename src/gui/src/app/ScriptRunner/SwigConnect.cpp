/**\file
 * Helpers to extract data from Python objects that are really SWIG wrappers
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include <Python.h>
#include "SwigConnect.h"

#include "SwigRuntime.h"

#include <iostream>

bool SwigConnect::SwigExtract (void*& data, PyObject* pyObj, const char* swigTypeStr)
{
  swig_type_info* type = SWIG_TypeQuery (swigTypeStr);
  if (!type)
  {
    // @@@ Proper error feedback (exception?)
    std::cerr << "couldn't find swig type info for: " << swigTypeStr << std::endl;
    return false;
  }
  
  if (SWIG_Python_ConvertPtr (pyObj, &data, type, 0) != 0)
  {
    // @@@ Proper error feedback (exception?)
    std::cerr << "SWIG_Python_ConvertPtr fail" << std::endl;
    return false;
  }
  if (!data)
  {
    // @@@ Proper error feedback (exception?)
    std::cerr << "SWIG_Python_ConvertPtr returned nullptr" << std::endl;
    return false;
  }
  
  return true;
}
