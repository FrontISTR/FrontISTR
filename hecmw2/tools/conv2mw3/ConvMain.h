/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   ConvMain.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef FA40620_1E0E_430b_BB4A_45051EC37AE8
#define FA40620_1E0E_430b_BB4A_45051EC37AE8

#include "ConvHeader.h"

class CConvMain{
private:
	CConvMain();
public:
	static CConvMain* Instance(){
		static CConvMain oConv;
		return &oConv;
	}
	virtual ~CConvMain();

private:
	CAssyModel moAssyModel;

	CFileReader moReader;
	CFileWriter moWriter;

public:
	void FileRead_FISTR4(string filename, bool opflag);// ファイル入力:FrontISTR v4 メッシュ
	void FileWrite_MW3(string filename);  // ファイル出力:MW3 メッシュ

	bool setupAssyModel();// AssyModelの出力準備(データ変換準備)
};
#endif //include_guard


