/* 
 * File:   LoggerType.h
 * Author: ktakeda
 *
 * Modify     2009/05/08
 * Created on 2009/05/01, 18:04
 */

#ifndef _LOGGERTYPE_H_ed5a0f37_3221_4434_8701_9ee348e993c4
#define	_LOGGERTYPE_H_ed5a0f37_3221_4434_8701_9ee348e993c4

namespace Utility{
// Logger Mode #
//
union LoggerMode{
enum{
    MWDebug,
    Debug,
    Info,
    Warn,
    Error,
    Invalid
};
};

// Logger Output Device #
//
union LoggerDevice{
enum{
    Disk,
    Display
};
};
}
#endif	/* _LOGGERTYPE_H */

