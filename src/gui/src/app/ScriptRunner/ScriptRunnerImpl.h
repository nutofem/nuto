#ifndef __SCRIPTRUNNERIMPL_H__
#define __SCRIPTRUNNERIMPL_H__

#include "ScriptRunner.h"
#include "ScriptRunnerThreaded.h"

class ScriptRunnerImpl : public nutogui::ScriptRunner
{
  ScriptRunnerThreaded threadedRunner;
public:
  /**\name nutogui::ScriptRunner implementation
   * @{ */
  void AddFeedbackCallback (const FeedbackCallbackPtr& callback);
  void StartScript (const wxString& source);
  /** @} */
};

#endif // __SCRIPTRUNNERIMPL_H__
