/**\file
 * Panel for controlling the scaling of displacement direction displaying.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "DisplacementDirectionSizePanel.h"

#include "uicommon/ColorFeedbackTextCtrl.h"
#include "uicommon/FloatingPointValidator.h"

#include "DirScaleDragThing.h"

#include <float.h>
#include <sstream>
#include <wx/artprov.h>

namespace nutogui
{
  DEFINE_EVENT_TYPE(EVENT_DIRECTION_SCALE_CHANGED)
  
  class DisplacementDirectionSizePanel::ScaleChangedEvent : public wxCommandEvent
  {
  public:
    ScaleChangedEvent () : wxCommandEvent (EVENT_DIRECTION_SCALE_CHANGED) {}
      
    wxEvent* Clone() const { return new ScaleChangedEvent (); }
  };
  
  //-------------------------------------------------------------------------
  
  enum
  {
    ID_ScaleInput = 1,
    ID_DragThing
  };
  
  BEGIN_EVENT_TABLE(DisplacementDirectionSizePanel, wxPanel)
    EVT_TEXT_ENTER (ID_ScaleInput, DisplacementDirectionSizePanel::OnScaleInputEnter)
    EVT_TEXT (ID_ScaleInput, DisplacementDirectionSizePanel::OnScaleInputChanged)
    
    EVT_COMMAND(ID_DragThing, EVENT_SCALEDRAG_BEGIN, DisplacementDirectionSizePanel::OnScaleDragBegin)
    EVT_SCALEDRAG_CHANGED(ID_DragThing, DisplacementDirectionSizePanel::OnScaleDragChanged)
    
    EVT_MOUSEWHEEL(DisplacementDirectionSizePanel::OnMouseWheel)
  END_EVENT_TABLE()
  
  DisplacementDirectionSizePanel::DisplacementDirectionSizePanel (wxWindow* parent, wxPoint position)
   : wxPanel (parent, wxID_ANY, position), locale (""), scaleValue (1),
     lastWheelTimeStamp (0), wheelDist (0)
  {
    wxColour bgColour (32, 32, 32);
    wxColour textColour (240, 240, 240);
    SetBackgroundColour (bgColour); // Somehow that results in a *different* color.
    
    wxBoxSizer* sizer = new wxBoxSizer (wxHORIZONTAL);
    
    wxStaticBitmap* decoBitmap = new wxStaticBitmap (this, wxID_ANY,
						     wxArtProvider::GetBitmap (wxART_MAKE_ART_ID (displace-scale-decoration)));
    sizer->Add (decoBitmap, wxSizerFlags().Center().Border (wxLEFT, 2));
    
    scaleInput = new uicommon::ColorFeedbackTextCtrl (this, ID_ScaleInput,
						      wxEmptyString, wxDefaultPosition, wxDefaultSize,
						      wxTE_RIGHT | wxBORDER_NONE | wxTE_PROCESS_ENTER);
    scaleInput->SetValidator (uicommon::FloatingPointValidator (scaleValue, locale));
    scaleInput->SetBackgroundColour (bgColour);
    scaleInput->SetForegroundColour (textColour);
    scaleInput->SetToolTip (wxT ("Enter a displacement scale value"));
    // Tweak size of scale input to allow containment of large numbers
    {
      wxScreenDC dc;
      dc.SetFont (scaleInput->GetFont ());
      
      std::stringstream outstream;
      outstream.imbue (locale);
      typedef std::num_put<char> NumPutter;
      std::use_facet<NumPutter> (locale).put (outstream, outstream, ' ', FLT_MAX);
      if (!outstream.fail())
      {
	wxString numberStr (outstream.str().c_str(), wxConvLocal);
	numberStr.Append (wxT(" ")); // for good measure
	wxSize numberStrExtent (dc.GetTextExtent (numberStr));
	wxSize scaleInputMinSize (scaleInput->GetMinSize ());
	scaleInputMinSize.x = std::max (scaleInputMinSize.x, numberStrExtent.x);
	scaleInput->SetMinSize (scaleInputMinSize);
      }
    }
    sizer->Add (scaleInput);
    
    DirScaleDragThing* dragThing = new DirScaleDragThing (this, ID_DragThing);
    dragThing->SetBackgroundColour (bgColour);
    dragThing->SetCursor (wxCursor (wxCURSOR_SIZENS));
    dragThing->SetToolTip (wxT ("Click and drag to change displacement scale"));
    sizer->Add (dragThing, wxSizerFlags().Expand());
    
    SetSizerAndFit (sizer);
    
    TransferDataToWindow();
  }

  void DisplacementDirectionSizePanel::SetDisplacementScale (float value)
  {
    scaleValue = value;
    TransferDataToWindow ();
  }

  void DisplacementDirectionSizePanel::OnScaleInputEnter (wxCommandEvent& event)
  {
    if (TransferDataFromWindow())
    {
      ScaleChangedEvent changedEvent;
      changedEvent.SetEventObject (this);
      wxPostEvent (this, changedEvent);
      
      scaleInput->ResetFeedbackColorFactor ();
    }
    else
    {
      scaleInput->SetFeedbackColorFactor (1);
    }
  }

  void DisplacementDirectionSizePanel::OnScaleInputChanged (wxCommandEvent& event)
  {
    if (TransferDataFromWindow())
    {
      ScaleChangedEvent changedEvent;
      changedEvent.SetEventObject (this);
      wxPostEvent (this, changedEvent);
      
      scaleInput->ResetFeedbackColorFactor();
    }
    else
    {
      // Dimmed error color - the user may not be finished yet...
      scaleInput->SetFeedbackColorFactor (0.5);
    }
  }

  void DisplacementDirectionSizePanel::OnScaleDragBegin (wxCommandEvent& event)
  {
    scaleDragStartScale = scaleValue;
  }
  
  void DisplacementDirectionSizePanel::OnScaleDragChanged (DirScaleDragThing::ChangedEvent& event)
  {
    SetRelativeScale (event.GetTotalScale());
  }

  void DisplacementDirectionSizePanel::SetRelativeScale (float scale)
  {
    double newScale = scaleDragStartScale * scale;
    if (newScale > FLT_MAX)
      newScale = FLT_MAX;
    else if (newScale < FLT_EPSILON)
      newScale = FLT_EPSILON;
    scaleValue = newScale;
    TransferDataToWindow();

    ScaleChangedEvent changedEvent;
    changedEvent.SetEventObject (this);
    wxPostEvent (this, changedEvent);
    
    scaleInput->ResetFeedbackColorFactor();
  }

  void DisplacementDirectionSizePanel::OnMouseWheel (wxMouseEvent& event)
  {
    /* Accumulate wheel "distance" as long as the user keeps scrolling
       (s/he presumably stops scrolling when not having scrolled for 1000 ms) */
    if ((event.GetTimestamp() - lastWheelTimeStamp) > 1000)
    {
      scaleDragStartScale = scaleValue;
      wheelDist = 0;
    }

    // Use the "distance" to compute a factor for the scale
    wheelDist += event.GetWheelRotation () / double (event.GetWheelDelta ());
    // similar formula as for dragging; but climbs faster (as wheel scrolling is "harder" than dragging)
    float scale = pow (1.2, wheelDist);
    
    SetRelativeScale (scale);
    
    lastWheelTimeStamp = event.GetTimestamp();
  }
  
} // namespace nutogui
