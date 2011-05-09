//
//	FileReader.h
//
//				2009.01.09
//				2008.12.09
//				k.Takeda
#ifndef FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C
#define FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C

#include "CommonFile.h"
#include "CommonStd.h"
#include "TypeDef.h"

#include "MeshFactory.h"

#include "FileBlockName.h"
#include "Logger.h"

namespace FileIO{
class CFileReader{
public:
    CFileReader();
    virtual ~CFileReader();
private:
    string msLine;//getLineSt()
protected:
    pmw::CMeshFactory *mpFactory;// Mesh setup
    Utility::CLogger *mpLogger;

    string& getLineSt(ifstream& ifs);// getline()-> char* -> string
    bool TagCheck(string& s_line, const char* ctag);

public:
    virtual void setFactory(pmw::CMeshFactory *pFactory);

//    virtual bool Read(ifstream &ifs, char* line)=0;
    virtual bool Read(ifstream &ifs, string& sline)=0;
};
}
#endif
