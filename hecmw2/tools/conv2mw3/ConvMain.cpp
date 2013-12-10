/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   ConvMain.cpp
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#include "ConvMain.h"

CConvMain::CConvMain()
{
}
CConvMain::~CConvMain()
{
}
//--
// ファイル入力
//--
void CConvMain::FileRead_FISTR4(string filename, bool opflag)
{
	moReader.ReadMesh_FISTR4(&moAssyModel, filename, opflag);
}
//--
// ファイル出力
//--
void CConvMain::FileWrite_MW3(string filename)
{
	moWriter.WriteMesh_MW3(&moAssyModel, filename);
}

//--
// AssyModelの出力準備(データ変換準備)
//--
bool CConvMain::setupAssyModel()
{
	return moAssyModel.setup();
}



