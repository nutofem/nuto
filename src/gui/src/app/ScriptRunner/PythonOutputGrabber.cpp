/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include <boost/python.hpp>
#include "common.h"
#include "PythonOutputGrabber.h"

using namespace boost::python;

void PythonOutputGrabber::BindToPython ()
{
  class_<PythonOutputGrabber> ("_PythonOutputGrabber",
			       no_init)
    .def("write", &PythonOutputGrabber::Write)
  ;
}

PythonOutputGrabber::PythonOutputGrabber (const OutputPendingNotifierPtr& callback)
 : bufferDirty (false), buffer (nullptr), bufferEnd (nullptr),
   bufferWritePos (nullptr), bufferUnread (nullptr),
   callback (callback)
{
}

PythonOutputGrabber::PythonOutputGrabber (const PythonOutputGrabber&)
{
  // Don't allow copies...
  assert (false);
}

size_t PythonOutputGrabber::BufferWriteable ()
{
  if (bufferUnread >= bufferWritePos)
    return bufferUnread - bufferWritePos;
  else
    return (bufferEnd - bufferWritePos) + (bufferUnread - buffer);
}

void PythonOutputGrabber::EnlargeBuffer (size_t extraChars)
{
  size_t bufferSize = bufferEnd - buffer;
  size_t newSize = std::min (bufferSize * 2, bufferSize + 1*1024*1024);
  if ((newSize-bufferSize) < extraChars)
  {
    newSize = bufferSize + ((extraChars + 1023) / 1024) * 1024;
  }
  size_t ofsUnread = bufferUnread-buffer;
  size_t ofsWritePos = bufferWritePos-buffer;
  buffer = (wchar_t*)realloc (buffer, newSize * sizeof (wchar_t));
  if (ofsUnread > ofsWritePos)
  {
    // Unread offset is after write pos, so move unread data to end of buffer
    size_t unreadSize = bufferSize - ofsUnread;
    size_t newUnreadOfs = newSize - unreadSize;
    memmove (buffer + newUnreadOfs, buffer + ofsUnread, unreadSize * sizeof (wchar_t));
    ofsUnread = newUnreadOfs;
  }
  bufferEnd = buffer + newSize;
  bufferWritePos = buffer + ofsWritePos;
  bufferUnread = buffer + ofsUnread;
}

void PythonOutputGrabber::Write (const std::string& str)
{
  bool wasDirty;
  {
    wxCriticalSectionLocker _lock (bufferLock);
    
    size_t wcharsNeeded = wxConvLibc.ToWChar (nullptr, 0, str.c_str(), str.length());
    wchar_t* convertDest = (wchar_t*)alloca (wcharsNeeded * sizeof (wchar_t));
    wxConvLibc.ToWChar (convertDest, wcharsNeeded, str.c_str(), str.length());
  #if wxABI_VERSION < 20900
    // Get rid of NUL
    wcharsNeeded -= 1;
  #endif
    
    size_t bufferWriteable = BufferWriteable();
    if (bufferWriteable < wcharsNeeded)
    {
      EnlargeBuffer (wcharsNeeded - bufferWriteable);
    }
    if (bufferWritePos < bufferUnread)
    {
      // There is a continuous chunk available at WritePos
      assert ((bufferUnread-bufferWritePos) >= long (wcharsNeeded));
      memcpy (bufferWritePos, convertDest, wcharsNeeded * sizeof (wchar_t));
      bufferWritePos += wcharsNeeded;
    }
    else
    {
      // May need to split buffer in two
      size_t chunk1 = std::min (wcharsNeeded, size_t (bufferEnd - bufferWritePos));
      memcpy (bufferWritePos, convertDest, chunk1 * sizeof (wchar_t));
      bufferWritePos += chunk1;
      if (bufferWritePos == bufferEnd)
	bufferWritePos = buffer;
      size_t chunk2 = wcharsNeeded - chunk1;
      if (chunk2 > 0)
      {
	assert ((bufferUnread-bufferWritePos) >= long (chunk2));
	memcpy (bufferWritePos, convertDest + chunk1, chunk2 * sizeof (wchar_t));
	bufferWritePos += chunk2;
      }
    }

    wasDirty = bufferDirty;
    bufferDirty = true;
  }
  if (!wasDirty)
    callback->OutputPending (this);
}

#if wxUSE_UNICODE
  #if ((wxABI_VERSION >= 20900) && defined (wxUSE_UNICODE_WCHAR)) \
    || ((wxABI_VERSION < 20900) && defined (wxUSE_WCHAR_T))
    #define WXSTRING_WCHAR
  #endif
#endif

#ifdef wxUSE_UNICODE_UTF8
  #define wxConvString	wxConvUTF8
#else
  #define wxConvString	wxConvLibc
#endif

void PythonOutputGrabber::GetPendingOutput (wxString& out)
{
  wxCriticalSectionLocker _lock (bufferLock);
  
  if (bufferWritePos >= bufferUnread)
  {
    size_t pendingSize = bufferWritePos - bufferUnread;
    if (pendingSize == 0)
    {
      out.Empty();
      return;
    }
    
  #ifdef WXSTRING_WCHAR
    wxStringBuffer outBuf (out, pendingSize + 1);
    memcpy (outBuf, bufferUnread, pendingSize*sizeof (wchar_t));
    outBuf[pendingSize] = 0;
  #else
    size_t charsNeeded = wxConvLibc.FromWChar (nullptr, 0, bufferUnread, pendingSize);
    char convertDest[charsNeeded];
    wxConvString.FromWChar (convertDest, charsNeeded, bufferUnread, pendingSize);
    wxStringBuffer outBuf (out, charsNeeded);
    memcpy (outBuf, convertDest, charsNeeded);
  #endif
    bufferUnread += pendingSize;
  }
  else
  {
    size_t pending1 = bufferEnd - bufferUnread;
    size_t pending2 = bufferWritePos - buffer;
    if (pending1 + pending2 == 0)
    {
      out.Empty();
      return;
    }
    
  #ifdef WXSTRING_WCHAR
    wxStringBuffer outBuf (out, pending1 + pending2 + 1);
    memcpy (outBuf, bufferUnread, pending1*sizeof (wchar_t));
    memcpy (outBuf, buffer, pending2*sizeof (wchar_t));
    outBuf[pending1 + pending2] = 0;
  #else
    wchar_t pendingCombined[pending1 + pending2];
    memcpy (pendingCombined, bufferUnread, pending1 * sizeof (wchar_t));
    memcpy (pendingCombined + pending1, buffer, pending2 * sizeof (wchar_t));
    size_t charsNeeded = wxConvLibc.FromWChar (nullptr, 0, pendingCombined, pending1 + pending2);
    char convertDest[charsNeeded];
    wxConvString.FromWChar (convertDest, charsNeeded, pendingCombined, pending1 + pending2);
    wxStringBuffer outBuf (out, charsNeeded);
    memcpy (outBuf, convertDest, charsNeeded);
  #endif
  
    bufferUnread = buffer + pending2;
  }
  
  bufferDirty = false;
}

