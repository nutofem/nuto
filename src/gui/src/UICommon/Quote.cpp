/**\file
 * Helpers for putting nice-looking quotes around strings.
 */
/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "uicommon/Quote.h"

#include <wx/intl.h>

namespace uicommon
{
#if wxABI_VERSION >= 20900
  typedef char wxTranslateChar;
#else
  typedef wxChar wxTranslateChar;
#endif
  
  static wxString QuoteString (const wxChar* str,
			       const wxTranslateChar* lquote, wchar_t lquoteDef,
			       const wxTranslateChar* rquote, wchar_t rquoteDef)
  {
    // TODO: Actually translate
    wxString lquoteTranslated (lquote);
    wxString rquoteTranslated (rquote);
    if ((lquoteTranslated == lquote) && (rquoteTranslated == rquote))
    {
      const wchar_t lqStr[2] = { lquoteDef, 0 };
      const wchar_t rqStr[2] = { rquoteDef, 0 };
      lquoteTranslated = lqStr;
      rquoteTranslated = rqStr;
    }
    
    wxString result;
    result.Alloc (lquoteTranslated.Length() + wxStrlen (str) + rquoteTranslated.Length() + 1);
    result = lquoteTranslated;
    result.Append (str);
    result.Append (rquoteTranslated);
    return result;
  }
  
  wxString Quote::Single (const wxChar* str)
  {
    return QuoteString (str,
			wxTRANSLATE("`"), 0x2018,
			wxTRANSLATE("'"), 0x2019);
  }
  
  wxString Quote::Double (const wxChar* str)
  {
    return QuoteString (str,
			wxTRANSLATE("``"), 0x201c,
			wxTRANSLATE("''"), 0x201d);
  }
  
} // namespace uicommon
