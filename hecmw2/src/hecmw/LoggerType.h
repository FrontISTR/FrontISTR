/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/LoggerType.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef _LOGGERTYPE_H_ed5a0f37_3221_4434_8701_9ee348e993c4
#define	_LOGGERTYPE_H_ed5a0f37_3221_4434_8701_9ee348e993c4
namespace Utility{
struct LoggerMode{
enum{
    MWDebug,
    Debug,
    Info,
    Warn,
    Error,
    Invalid
};
};
struct LoggerDevice{
enum{
    Disk,
    Display
};
};
}
#endif	/* _LOGGERTYPE_H */
