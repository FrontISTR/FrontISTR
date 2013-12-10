/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   FileReader.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef A7595_5F38_4409_8F3B_A96BA007F4E5
#define A7595_5F38_4409_8F3B_A96BA007F4E5

#include "FileBase.h"
#include "AssyModel.h"

class CFileReader:public CFileBase{
public:
	CFileReader();
	virtual ~CFileReader();

protected:
	size_t mnLineNum;    //読み込み行数
	string msHeaderTitle;//!HEADER Title:ファイルの内容説明文

	bool mbOpFlag;//オプション・フラグ

protected:
	size_t getElementNNum(size_t nType);//要素タイプ別の節点数

	CMesh* getgenMesh(CAssyModel* pAssyModel, string sPartsName, string sTagName);
	CGroupPair* getgenMPCPair(CAssyModel* pAssyModel, string sPairName, string sTagName);
	CGroupPair* getgenConPair(CAssyModel* pAssyModel, string sPairName, string sTagName);
	CNgroup* getgenNgrp(CMesh* pMesh, string sNgrpName, string sTagName);
	CSgroup* getgenSgrp(CMesh* pMesh, string sSgrpName, string sTagName);
	CEgroup* getgenEgrp(CMesh* pMesh, string sEgrpName, string sTagName);
	CLgroup* getgenLgrp(CMesh* pMesh, string sLgrpName, string sTagName);
	
	CElement* getgenElement(CMesh* pMesh, size_t nID, size_t nType);

	bool ReadHeader(ifstream& ifs, string& sLine);
	bool ReadNode(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);
	bool ReadElement(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);
	bool ReadEgrp(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);
	bool ReadSgrp(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);
	bool ReadNgrp(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);
	bool ReadGrpPair(CAssyModel *pAssyModel, ifstream& ifs, string& sLine);//Assembly_Pair &Contact_Pair

public:
	void ReadMesh_FISTR4(CAssyModel *pAssyModel, string filename, bool opflag);
};
#endif //include_guard





