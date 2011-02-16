/**\file
 * Interface for dealing with script errors
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_SCRIPTERROR_H__
#define __NUTOGUI_SCRIPTERROR_H__

#include "ScriptRunner.h"

#include <boost/shared_ptr.hpp>

namespace nutogui
{
  struct ScriptError
  {
    /// Clear current error location (i.e. remove error indicators)
    virtual void ClearErrorLocation () = 0;
    /**
     * Set current error location from traceback.
     * \return Whether the error could be located (i.e. it happened in a known script).
     */
    virtual bool SetErrorLocation (const ScriptRunner::TracebackPtr& traceback) = 0;
    /// Display error location (move caret to it etc.)
    virtual void DisplayErrorLocation () = 0;
  };
  typedef boost::shared_ptr<ScriptError> ScriptErrorPtr;
} // namespace nutogui

#endif // __NUTOGUI_SCRIPTERROR_H__
