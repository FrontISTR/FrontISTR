/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Sgroup.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef DC03BBA7_5463_4bc3_889E_E8202FDE0E98
#define DC03BBA7_5463_4bc3_889E_E8202FDE0E98

#include "Group.h"

class CSgroup:public CGroup{
public:
	CSgroup();
	virtual ~CSgroup();

private:
	vector<pair<size_t,size_t> > mvElemFace;//入力用:要素番号-局所面番号ペア(FrontISTR v4)


	vector<pair<size_t,size_t> > mvElemFaceMW3;//要素番号-局所面番号ペア(MW3)
	vector<CNode*> mvNode;       //Sgroupに所属するNode*
	map<size_t,size_t> mmNID2NIX;//Node_id⇒index

	vector<size_t> mvFaceNodeNum;//面の節点数
	vector<string> mvFaceType;   //面のタイプ
	vector<size_t> mvMW3FaceNum; //面番号(MW3)
	vector<vector<size_t> > mvvBNodeNum;//面構成のBNode番号(index番号)
	vector<vector<size_t> > mvvFaceNodeID;//面構成のNodeID

public:
	void addElemFaceID(size_t nElemID, size_t nFaceN);

	//--
	// DATA出力準備
	//--
	bool setup(map<size_t,CNode*> mNode, map<size_t,CElement*> mElem);

	//--
	// 出力
	//--
	size_t getNumOfNode();
	CNode* getNode(size_t index);

	size_t getNumOfFace();
	size_t getElemID(size_t index);
	size_t getFaceN(size_t index);//面番号 入力値のまま.

	size_t getFaceNodeNum(size_t index);
	string getFaceType(size_t index);
	size_t getConvMW3FaceN(size_t index);//FrontISTR⇒MW3 面番号
	vector<size_t> getBNodeN_Face(size_t index);//面を構成するBoundaryNode番号(mvNodeのindex番号)

	vector<size_t> getFaceNodeID(size_t index);//生データ(入力値)

	map<size_t,size_t> getNID2NIX(){ return mmNID2NIX;}
	
};
#endif //include_guard


