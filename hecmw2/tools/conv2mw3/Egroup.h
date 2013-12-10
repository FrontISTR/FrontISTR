/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Egroup.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef F48A2364_80F8_41c8_912E_6C1FD23F9982
#define F48A2364_80F8_41c8_912E_6C1FD23F9982

#include "Element.h"
#include <vector>

#include "Group.h"

class CEgroup:public CGroup{
public:
	CEgroup();
	virtual ~CEgroup();

private:
	vector<CElement*> mvTempElem;
	vector<CElement*> mvElement;
	map<size_t,size_t> mmEID2IX;
	vector<size_t> mvEID;

	vector<size_t> mvNID;//EgrpëSëÃÇÃNodeIDî‘çÜîzóÒ
	map<size_t,size_t> mmNID2IX;//NodeIDÅÀmvNIDÇÃIndex

public:
	void addElement(CElement* pElem);

	size_t getNumOfElement(){return mvElement.size();}
	CElement* getElement(size_t index){ return mvElement[index];}

	//--
	// DATAèoóÕèÄîı
	//--
	bool setup();

	size_t getNumOfNode();
	size_t getNodeID(size_t index);
	size_t getBNodeIndex(size_t nid);
};

#endif //include_guard



