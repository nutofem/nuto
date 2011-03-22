/**\file
 * Control to change the scale of displacement direction display by dragging the mouse.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"

#include "DirScaleDragThing.h"

#include <wx/artprov.h>
#include <wx/dcbuffer.h>

namespace nutogui
{
  DEFINE_EVENT_TYPE(EVENT_SCALEDRAG_BEGIN)
  DEFINE_EVENT_TYPE(EVENT_SCALEDRAG_CHANGED)
  DEFINE_EVENT_TYPE(EVENT_SCALEDRAG_END)
  
  enum
  {
    TIMER_HotFade = 1
  };
  
  BEGIN_EVENT_TABLE (DirScaleDragThing, wxControl)
    EVT_PAINT (DirScaleDragThing::OnPaint)
    
    EVT_TIMER (TIMER_HotFade, DirScaleDragThing::OnHotFadeTimer)
    
    EVT_LEFT_DOWN(DirScaleDragThing::OnMouseLeftDown)
    EVT_LEFT_UP(DirScaleDragThing::OnMouseLeftUp)
    EVT_MOTION(DirScaleDragThing::OnMouseMove)
    EVT_MOUSE_CAPTURE_LOST(DirScaleDragThing::OnMouseCaptureLost)
    EVT_ENTER_WINDOW(DirScaleDragThing::OnMouseEnter)
    EVT_LEAVE_WINDOW(DirScaleDragThing::OnMouseLeave)
  END_EVENT_TABLE()

  static const int hotFadeOutTime = 500;
  
  DirScaleDragThing::DirScaleDragThing (wxWindow* parent, wxWindowID id,
					wxPoint position, wxSize size)
   : wxControl (parent, id, position, size, wxBORDER_NONE),
     isHot (false), hotFadeTimer (this, TIMER_HotFade)
  {
    SetBackgroundStyle (wxBG_STYLE_CUSTOM);
    
    image = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID (displace-scale-drag-thing));
    hotImage = wxArtProvider::GetBitmap (wxART_MAKE_ART_ID (displace-scale-drag-thing-hot));
    
    SetMinSize (wxSize (image.GetWidth(), image.GetHeight()));
  }

  void DirScaleDragThing::OnPaint (wxPaintEvent& event)
  {
    wxAutoBufferedPaintDC dc (this);
    float fadeFactor = isHot ? 1 : (1 - (hotFadeStopWatch.Time() / float (hotFadeOutTime)));
    PaintImage (dc, std::max (fadeFactor, 0.0f));
  }

  void DirScaleDragThing::OnMouseLeftDown (wxMouseEvent& event)
  {
    dragStart = event.GetPosition();
    CaptureMouse ();
    
    wxCommandEvent dragStartEvent (EVENT_SCALEDRAG_BEGIN, GetId());
    dragStartEvent.SetEventObject (this);
    GetEventHandler()->ProcessEvent (dragStartEvent);
  }

  void DirScaleDragThing::OnMouseLeftUp (wxMouseEvent& event)
  {
    if (HasCapture ())
    {
      ReleaseMouse ();
      DragEnd ();
      SetHot (false);
    }
  }

  void DirScaleDragThing::OnMouseMove (wxMouseEvent& event)
  {
    if (HasCapture ())
    {
      wxPoint newPos (event.GetPosition());
      
      // Distance to a line through the drag start, goint from bottom left to top right
      int distToStart = (newPos.x + newPos.y - dragStart.x - dragStart.y);
      // Exact distance would be * M_SQRT1_2 ... we don't need that
      
      float scale = pow (1.05, distToStart);
      ChangedEvent changeEvent (GetId(), scale);
      changeEvent.SetEventObject (this);
      GetEventHandler()->ProcessEvent (changeEvent);
    }
  }

  void DirScaleDragThing::OnMouseCaptureLost (wxMouseCaptureLostEvent& event)
  {
    DragEnd ();
    SetHot (false);
  }
    
  void DirScaleDragThing::OnMouseEnter (wxMouseEvent& event)
  {
    SetHot (true);
  }

  void DirScaleDragThing::OnMouseLeave (wxMouseEvent& event)
  {
    if (!HasCapture())
    {
      SetHot (false);
    }
  }
  
  void DirScaleDragThing::DragEnd ()
  {
    wxCommandEvent dragEndEvent (EVENT_SCALEDRAG_END, GetId());
    dragEndEvent.SetEventObject (this);
    GetEventHandler()->ProcessEvent (dragEndEvent);
  }

  void DirScaleDragThing::OnHotFadeTimer (wxTimerEvent& event)
  {
    float fadeFactor = 1 - (hotFadeStopWatch.Time() / float (hotFadeOutTime));
    if (fadeFactor <= 0)
    {
      hotFadeTimer.Stop();
    }
    
    wxClientDC clientDC (this);
    if (IsDoubleBuffered())
    {
      PaintImage (clientDC, fadeFactor);
    }
    else
    {
      wxBufferedDC dc (&clientDC, GetSize());
      PaintImage (dc, fadeFactor);
    }
  }
    
  void DirScaleDragThing::PaintImage (wxDC& dc, float fade)
  {
    dc.SetPen (*wxTRANSPARENT_PEN);
    dc.SetBrush (wxBrush (GetBackgroundColour ()));
    dc.DrawRectangle (0, 0, GetSize().GetWidth(), GetSize().GetHeight());
    
    wxCoord x = (GetSize().GetWidth() - image.GetWidth()) / 2;
    wxCoord y = (GetSize().GetHeight() - image.GetHeight()) / 2;
    if (fade > 0)
    {
      if (fade < 1)
      {
	DrawBitmapsMixed (dc, x, y, image, hotImage, fade);
      }
      else
	dc.DrawBitmap (hotImage, x, y);
    }
    else
      dc.DrawBitmap (image, x, y);
  }

  void DirScaleDragThing::SetHot (bool hotState)
  {
    if (hotState)
    {
      isHot = true;
      hotFadeTimer.Stop ();
      Refresh (false);
    }
    else
    {
      if (!isHot) return;
      isHot = false;
      hotFadeStopWatch.Start();
      hotFadeTimer.Start (50);
      Refresh (false);
    }
  }

  void DirScaleDragThing::DrawBitmapsMixed (wxDC& dc, wxCoord x, wxCoord y,
					    const wxBitmap& image1, const wxBitmap& image2,
					    float alpha)
  {
    /* @@@ This is all pretty roundabout, but there doesn't seem to be a more
       direct way, at least in WX 2.8 */
    assert ((image1.GetWidth() == image2.GetWidth()) && (image1.GetHeight() == image2.GetHeight()));
    
    wxImage img1 (image1.ConvertToImage());
    if (!img1.HasAlpha()) img1.InitAlpha();
    wxImage img2 (image2.ConvertToImage());
    if (!img2.HasAlpha()) img2.InitAlpha();
    
    wxImage destImage (image1.GetWidth(), image1.GetHeight(), false);
    destImage.InitAlpha();
    
    float factor1 = 1-alpha;
    float factor2 = alpha;
    const size_t numPixels = image1.GetWidth() * image1.GetHeight();
    {
      const unsigned char* rgb1 = img1.GetData();
      const unsigned char* a1 = img1.GetAlpha();
      const unsigned char* rgb2 = img2.GetData();
      const unsigned char* a2 = img2.GetAlpha();
      unsigned char* rgbDst = destImage.GetData();
      unsigned char* aDst = destImage.GetAlpha();
      for (size_t p = 0; p < numPixels; p++)
      {
	for (int c = 0; c < 3; c++)
	{
	  *rgbDst = (*rgb1)*factor1 + (*rgb2)*factor2;
	  rgbDst++;
	  rgb1++;
	  rgb2++;
	}
	*aDst = (*a1)*factor1 + (*a2)*factor2;
	aDst++;
	a1++;
	a2++;
      }
    }
    
    wxBitmap bmp (destImage);
    dc.DrawBitmap (bmp, x, y);
  }
  
  //-------------------------------------------------------------------------
  
  DirScaleDragThing::ChangedEvent::ChangedEvent (wxWindowID id, float totalScale)
   : wxCommandEvent (EVENT_SCALEDRAG_CHANGED, id), totalScale (totalScale) {}
  
  wxEvent* DirScaleDragThing::ChangedEvent::Clone() const
  {
    return new ChangedEvent (*this);
  }

} // namespace nutogui
