/*
 ----------------------------------------------------------
|
| Software Name :HEC middleware Ver. 3.0beta
|
|   FileReader.h
|
|                     Written by T.Takeda,    2010/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C
#define FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C
#include "CommonFile.h"
#include "CommonStd.h"
#include "TypeDef.h"
#include "MeshFactory.h"
#include "FileBlockName.h"
#include "Logger.h"
#include "boost/lexical_cast.hpp"
#include "ElementType.h"
namespace FileIO{
class CFileReader{
public:
    CFileReader();
    virtual ~CFileReader();
private:
    string msLine;
protected:
    pmw::CMeshFactory *mpFactory;
    Utility::CLogger *mpLogger;
    string& getLineSt(ifstream& ifs);
    bool TagCheck(string& s_line, const char* ctag);
    uint IntElemType(string& sElemType);
public:
    virtual void setFactory(pmw::CMeshFactory *pFactory);
    virtual bool Read(ifstream &ifs, string& sline)=0;
};
}
#endif
