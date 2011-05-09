//
//  ContactMesh.h
//
//
//				2009.01.08
//				2009.01.08
//				k.Takeda
#ifndef CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F
#define CONTACT_MESH_HH_958EE36F_3904_4b1e_91B5_F6308F5AAC9F

#include "CommonStd.h"
#include "ContactElement.h"

namespace pmw{
class CContactMesh{
public:
	CContactMesh();
	virtual ~CContactMesh();

protected:
	//
	uint mnContactID;

	// Contact Element (Face) 
	vector<CContactElement*> mvContactElement;

public:
	void setID(const uint& id){ mnContactID = id;}
	uint& getID(){ return mnContactID;}

	uint getNumOfContactElem(){ return mvContactElement.size();}

	// Face.push_back
	void addContactElement(CContactElement* pContactElem){ mvContactElement.push_back(pContactElem);}

	//
	void resizeContactElem(const uint& size){ mvContactElement.resize(size);}
	void setContactElem(CContactElement* pContactElem, const uint& i){ mvContactElement[i] = pContactElem;}
};
}
#endif
