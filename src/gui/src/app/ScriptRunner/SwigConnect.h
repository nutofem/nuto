/**\file
 * Helpers to extract data from Python objects that are really SWIG wrappers
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __SWIGCONNECT_H__
#define __SWIGCONNECT_H__

struct SwigConnect
{
public:
  static bool SwigExtract (void*& data, PyObject* pyObj, const char* swigTypeStr);
  
  template<typename T>
  static bool SwigExtract (T*& data, PyObject* pyObj, const char* swigTypeStr)
  {
    void* p;
    if (!SwigExtract (p, pyObj, swigTypeStr)) return false;
    data = reinterpret_cast<T*> (p);
    return true;
  }
};

#endif // __SWIGCONNECT_H__
