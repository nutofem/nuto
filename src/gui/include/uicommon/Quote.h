/**\file
 * Helpers for putting nice-looking quotes around strings.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_QUOTE_H__
#define __UICOMMON_QUOTE_H__

#include "export.h"

#include <wx/string.h>

namespace uicommon
{
  /// Helpers for putting nice-looking quotes around strings
  struct UICOMMON_EXPORTED Quote
  {
    //@{
    /// Put single quotes around a string
    static wxString Single (const wxChar* str);
    static wxString Single (const wxString& str)
    {
      return Single ((const wxChar*)str);
    }
    //@}
    
    //@{
    /// Put single quotes around a string
    static wxString Double (const wxChar* str);
    static wxString Double (const wxString& str)
    {
      return Double ((const wxChar*)str);
    }
    //@}
  };
} // namespace uicommon

#endif // __UICOMMON_QUOTE_H__
