/*
 ----------------------------------------------------------
|
| Software Name :conv2mw3 Ver 0.1 beta
|
|   Node.h
|
|                     Written by T.Takeda,    2012/06/01
|                                Y.Sato,      2012/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
#ifndef B0F859_DE70_4dee_88CE_7F57C05DF4DB
#define B0F859_DE70_4dee_88CE_7F57C05DF4DB

#include <cstdlib>

class CNode{
public:
	CNode(void);
	virtual ~CNode(void);

private:
	size_t mID;
	double mdCoord[3];

public:
	void setID(size_t id);
	void setCoord(const double& x, const double& y, const double& z);

	size_t getID(){return mID;}

	double& getX(){return mdCoord[0];}
	double& getY(){return mdCoord[1];}
	double& getZ(){return mdCoord[2];}

};

#endif //include guard 
