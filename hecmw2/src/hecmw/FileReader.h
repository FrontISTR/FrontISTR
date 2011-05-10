//
//	FileReader.h
//
//				2009.01.09
//				2008.12.09
//				k.Takeda
#ifndef FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C
#define FILE_READER_HH_D67C7697_52EC_426c_A2D3_FB67B60AF90C

#include "CommonFile.h"
#include "TypeDef.h"

#include "MeshFactory.h"

#include "FileBlockName.h"
#include "Logger.h"

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include "ElementType.h"
#include "ElementProperty.h"

namespace FileIO{
class CFileReader{
public:
    CFileReader();
    virtual ~CFileReader();
private:
    string msLine;//getLineSt()
protected:
    pmw::CMeshFactory *mpFactory;// Mesh setup
    Utility::CLogger  *mpLogger;

    string& getLineSt(ifstream& ifs);// getline()-> char* -> string
    string& getLine(ifstream& ifs);// \r, \n, ',' => ' '
    bool TagCheck(string& s_line, const char* ctag);

    uint IntElemType(string& sElemType);
    uint IntBndType(string& sBndType);//文字列のDirichlet.or.Neumannをuintへ変更.

    void Split(const string& s, char c, vstring& v);// s:対象文字列、c:区切り文字 , v:分割されたトークン
public:
    virtual void setFactory(pmw::CMeshFactory *pFactory);

//    virtual bool Read(ifstream &ifs, char* line)=0;
    virtual bool Read(ifstream &ifs, string& sline)=0;
};
}
#endif
