/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Group.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef C538270_9FA4_4bcc_B1C7_40D803AF26FC
#define C538270_9FA4_4bcc_B1C7_40D803AF26FC

#include <utility>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
using namespace std;

#include "Node.h"
#include "Element.h"

#include "FrontISTR_Type.h"

class CGroup{
public:
	CGroup();
	virtual ~CGroup();

protected:
	size_t mID;
	string mGroupName;
	string mPartsName;

	size_t mnBndType;//ã´äEèåèéÌóﬁ BndType::NotUse
public:
	void setID(size_t id);
	size_t getID();

	void setGroupName(string name);
	void setPartsName(string name);

	string getGroupName(){ return mGroupName;}
	string getPartsName(){ return mPartsName;}

	void   setBndType(size_t nType);
	size_t getBndType();
};
#endif //include_guard









