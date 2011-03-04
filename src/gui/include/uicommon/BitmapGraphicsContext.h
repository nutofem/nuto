/**\file
 * Helper class to get a wxGraphicsContext for a wxBitmap
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __UICOMMON_BITMAPGRAPHICSCONTEXT_H__
#define __UICOMMON_BITMAPGRAPHICSCONTEXT_H__

#include "export.h"

#include <wx/dcmemory.h>

class wxBitmap;
class wxGraphicsContext;

namespace uicommon
{
#if wxUSE_GRAPHICS_CONTEXT
  typedef wxGraphicsContext GraphicsContext;
#else
  typedef wxDC GraphicsContext;
#endif
  
  /**
   * Helper class to get a wxGraphicsContext for a wxBitmap.
   */
  class UICOMMON_EXPORTED BitmapGraphicsContext
  {
    /**
     * Flag to hack around wxGraphicsContext for a wxMemoryDC with a wxBitmap
     * not copying the alpha data when updating the wxBitmap.
     */
    bool alphaHack;
    
    GraphicsContext* theGC;
    wxBitmap& bmp;
    wxMemoryDC dc;
  public:
    enum ClearMode { clear, dontClear };
    /**
     * Set up helper.
     * \param bmp Bitmap to draw on.
     * \param clearMode Whether to clear the bitmap (clear) or leave it as it is (dontClear)
     */
    BitmapGraphicsContext (wxBitmap& bmp, ClearMode clearMode = clear);
    ~BitmapGraphicsContext();
    
    GraphicsContext* operator->() { return theGC; }
  };
} // namespace uicommon

#endif // __UICOMMON_BITMAPGRAPHICSCONTEXT_H__
