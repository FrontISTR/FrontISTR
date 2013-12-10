/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   FileWriter.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef E22CB_31A8_46bd_BB46_8E538E03042C
#define E22CB_31A8_46bd_BB46_8E538E03042C

#include "FileBase.h"
#include "AssyModel.h"

#include "FrontISTR_Type.h"

class CFileWriter:public CFileBase{
public:
	CFileWriter();
	virtual ~CFileWriter();

private:
	string msDirichlet, msNeumann;

	//--
	// ブロック別・出力
	//--
	bool WriteHeader(ofstream& ofs);//-- Header
	bool WriteAssyModel(CAssyModel *pAssyModel, ofstream& ofs);//---AssyModel
	bool WriteNode(CAssyModel *pAssyModel, ofstream& ofs);     //---Node
	bool WriteElement(CAssyModel *pAssyModel, ofstream& ofs);  //---Element
	bool WriteBndNodeMesh(CAssyModel *pAssyModel, ofstream& ofs);//-BoundaryNodeMesh
	bool WriteBndFaceMesh(CAssyModel *pAssyModel, ofstream& ofs);//-BoundaryFaceMesh
	bool WriteBndEdgeMesh(CAssyModel *pAssyModel, ofstream& ofs);//-BoundaryEdgeMesh
	bool WriteBndVolMesh(CAssyModel *pAssyModel, ofstream& ofs); //-BoundaryVolumeMesh
	bool WriteElemGroup(CAssyModel *pAssyModel, ofstream& ofs);//---ElementGroup ← BoundaryVolumeMeshの読み替え
	//
	// # 通信界面はシングルなので無い : CommMesh2
	//
	bool WriteContactMesh(CAssyModel *pAssyModel, ofstream& ofs);//-ContactMesh

public:
	void WriteMesh_MW3(CAssyModel *pAssyModel, string filename);

};
#endif //include_guard




