/*
 * Written 2010 by Frank Richter <frank.richter@gmail.com>
 * For Bauhaus University Weimar, Institute of Structural Mechanics
 *
 */

#include "common.h"
#include "platform/ThreadRPC.h"

#include <wx/app.h>

namespace platform
{
  
class ThreadRPC::MethodPendingEvent : public wxEvent
{
  MainThreadRunnerBasePtr runner;
public:
  MethodPendingEvent (MainThreadRunnerBasePtr runner)
    : wxEvent (wxID_ANY, trpcEVT_METHOD_PENDING), runner (runner) {}
  
  wxEvent* Clone() const { return new MethodPendingEvent (runner); }
  
  MainThreadRunnerBasePtr GetRunner() const { return runner; }
};

//---------------------------------------------------------------------------
  
void ThreadRPC::ConditionSignalMain::Signal()
{
  MethodPendingEvent event (mtr->GetShared());
  wxPostEvent (wxTheApp, event);
}

//---------------------------------------------------------------------------
  
DEFINE_EVENT_TYPE(trpcEVT_METHOD_PENDING)

class MainThreadEventHandler : public wxEvtHandler
{
public:
  void OnMethodPending(ThreadRPC::MethodPendingEvent& event)
  {
    while (event.GetRunner()->Dispatch());
  }
};
typedef void (wxEvtHandler::*MethodPendingEventFunction)(ThreadRPC::MethodPendingEvent&);

MainThreadEventHandler mainThreadEventHandler;

static bool handlerInstalled = false;

void ThreadRPC::InstallMainThreadEventHandler()
{
  if (!handlerInstalled)
  {
    handlerInstalled = true;
    
    wxTheApp->Connect (trpcEVT_METHOD_PENDING,
		       (wxObjectEventFunction)(wxEventFunction)static_cast<MethodPendingEventFunction> (&MainThreadEventHandler::OnMethodPending),
		       nullptr,
		       &mainThreadEventHandler);
  }
}

} // namespace platform
