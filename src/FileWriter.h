/* 
 * File:   FileWriter.h
 * Author: ktakeda
 *
 * Created on 2009/07/23, 11:55
 */

#ifndef _FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e
#define	_FILEWRITER_H_6e7d31b0_0854_4c88_af18_efc72922f25e

#include "CommonFile.h"

#include "GMGModel.h"


namespace FileIO{
class CFileWriter{
public:
    CFileWriter();
    virtual ~CFileWriter();

protected:
    pmw::CGMGModel *mpGMGModel;

public:
    virtual void Write(ofstream& ofs, const uint& mgLevel)=0;
};
}
#endif	/* _FILEWRITER_H */

