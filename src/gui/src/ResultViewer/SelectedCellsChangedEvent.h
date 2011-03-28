/**\file
 * Event sent when the set of selected cells changes.
 */
/*
 * Written 2011 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#ifndef __NUTOGUI_RESULTVIEWER_SELECTEDCELLSCHANGEDEVENT_H__
#define __NUTOGUI_RESULTVIEWER_SELECTEDCELLSCHANGEDEVENT_H__

#include "ResultViewerImpl.h"

namespace nutogui
{
  BEGIN_DECLARE_EVENT_TYPES()
    DECLARE_EVENT_TYPE(EVENT_SELECTEDCELLS_CHANGED, 0)
  END_DECLARE_EVENT_TYPES()
  
  #define EVT_SELECTEDCELLS_CHANGED(fn)	\
    EVT_COMMAND(wxID_ANY, EVENT_SELECTEDCELLS_CHANGED, fn)

  class ResultViewerImpl::SelectedCellsChangedEvent : public wxEvent
  {
  public:
    SelectedCellsChangedEvent () : wxEvent (0, EVENT_SELECTEDCELLS_CHANGED) {}
      
    wxEvent* Clone() const { return new SelectedCellsChangedEvent (); }
  };
  
} // namespace nutogui

#endif // __NUTOGUI_RESULTVIEWER_SELECTEDCELLSCHANGEDEVENT_H__
