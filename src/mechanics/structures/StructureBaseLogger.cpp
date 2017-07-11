// $Id:$
#include "mechanics/structures/StructureBase.h"

//*************************************************
//************    Logger routines    ***************
//*************************************************
//! @brief opens a logger file for the output to a log file
//! @param rFileName file name
void NuTo::StructureBase::LoggerOpenFile(std::string rFileName)
{
    mLogger.OpenFile(rFileName);
}

//! @brief set the logger to be quiet (output only to file, if set)
//! @param rQuiet (true for quiet logger, false for output to standard output)
void NuTo::StructureBase::LoggerSetQuiet(bool rQuiet)
{
    mLogger.SetQuiet(rQuiet);
}
