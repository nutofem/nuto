#include "common.h"
#include "ScriptRunnerImpl.h"

#include <wx/wx.h>

void ScriptRunnerImpl::AddFeedbackCallback (const FeedbackCallbackPtr& callback)
{
  threadedRunner.AddFeedbackCallback (callback);
}

void ScriptRunnerImpl::StartScript (const wxString& source)
{
  threadedRunner.StartScript (source);
}
