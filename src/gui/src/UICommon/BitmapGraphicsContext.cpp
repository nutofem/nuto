/**\file
 * Helper class to get a wxGraphicsContext for a wxBitmap
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "uicommon/BitmapGraphicsContext.h"

#include <wx/graphics.h>
#include <wx/rawbmp.h>

namespace uicommon
{
  BitmapGraphicsContext::BitmapGraphicsContext (wxBitmap& bmp, ClearMode clearMode)
    : theGC (nullptr), bmp (bmp)
  {
    alphaHack = false;
    if (clearMode == clear)
    {
    #if defined(__WXGTK__)
      alphaHack = true;
    #endif
    
      // Clear bitmap data.
      wxAlphaPixelData bmpPixels (bmp);
      wxAlphaPixelData::Iterator p (bmpPixels);
      for (size_t numPix = bmpPixels.GetWidth() * bmpPixels.GetHeight(); numPix-- > 0; )
      {
	p.Red() = alphaHack ? 1 : 0;
	p.Green() = alphaHack ? 1 : 0;
	p.Blue() = alphaHack ? 1 : 0;
      #if wxUSE_GRAPHICS_CONTEXT
	p.Alpha() = 0;
      #else
	/* WX 2.8: setting the alpha to 0 and to 255 later doesn't always work...
	   whyever, whatever, invert the alpha hack logic in that case. */
	p.Alpha() = 255;
      #endif
	++p;
      }
    }
    
    dc.SelectObject (bmp);
  #if wxUSE_GRAPHICS_CONTEXT
    theGC = wxGraphicsContext::Create (dc);
  #else
    theGC = &dc;
  #endif
  }

  BitmapGraphicsContext::~BitmapGraphicsContext()
  {
  #if wxUSE_GRAPHICS_CONTEXT
    delete theGC;
  #endif
    dc.SelectObject (wxNullBitmap);
    
    if (alphaHack)
    {
      // Horrible, horrible - create *some* alpha by using a mask color
      wxAlphaPixelData bmpPixels (bmp);
      wxAlphaPixelData::Iterator p (bmpPixels);
      for (int y = 0; y < bmpPixels.GetHeight(); y++)
      {
	for (int x = 0; x < bmpPixels.GetWidth(); x++)
	{
	#if wxUSE_GRAPHICS_CONTEXT
	  if ((p.Red() != 1) || (p.Green() != 1) || (p.Blue() != 1))
	    p.Alpha() = 255;
	#else
	  // See above
	  if ((p.Red() == 1) && (p.Green() == 1) && (p.Blue() == 1))
	    p.Alpha() = 0;
	#endif
	++p;
	}
      }
    }
  }
  
} // namespace uicommon
