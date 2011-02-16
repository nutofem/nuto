/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __PYTHONOUTPUTGRABBER_H__
#define __PYTHONOUTPUTGRABBER_H__

#include <wx/string.h>
#include <wx/thread.h>

#include <boost/shared_ptr.hpp>
#include <string>

class PythonOutputGrabber
{
  wxCriticalSection bufferLock;
  bool bufferDirty;
  // Ring buffer for grabbed output
  wchar_t* buffer;
  wchar_t* bufferEnd;
  wchar_t* bufferWritePos;
  wchar_t* bufferUnread;
  
  // Return the number of chars that can be written to the buffer
  inline size_t BufferWriteable ();
  // Enlarge buffer so it holds (at least) 'extraChars' additional characters
  void EnlargeBuffer (size_t extraChars);
public:
  /// Register with Boost.Python
  static void BindToPython ();
  
  struct OutputPendingNotifier
  {
    virtual ~OutputPendingNotifier () {}
    
    virtual void OutputPending (PythonOutputGrabber* source) = 0;
  };
  typedef boost::shared_ptr<OutputPendingNotifier> OutputPendingNotifierPtr;
  
  PythonOutputGrabber (const OutputPendingNotifierPtr& callback);
  // Need to provide copy ctor so boost::python bindings build
  PythonOutputGrabber (const PythonOutputGrabber& other);
  
  /// Write something to buffer. Called from Python
  void Write (const std::string& str);
  
  /// Get all pending output, copy to \a out
  void GetPendingOutput (wxString& out);
private:
  OutputPendingNotifierPtr callback;
};

#endif // __PYTHONOUTPUTGRABBER_H__
