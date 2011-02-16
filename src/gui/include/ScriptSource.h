/**\file
 * Interface for script sources
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_SCRIPTSOURCE_H__
#define __NUTOGUI_SCRIPTSOURCE_H__

#include <wx/string.h>

#include <boost/shared_ptr.hpp>

namespace nutogui
{
  struct ScriptSource
  {
    virtual wxString GetSourceName() const = 0;
    virtual wxString GetSourceString() const = 0;
  };
  typedef boost::shared_ptr<ScriptSource> ScriptSourcePtr;
} // namespace nutogui

#endif // __NUTOGUI_SCRIPTSOURCE_H__
