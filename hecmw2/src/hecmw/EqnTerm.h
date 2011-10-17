/*
 ----------------------------------------------------------
|
| Software Name :HEC-MW Ver 4.0beta
|
|   ../src/EqnTerm.h
|
|                     Written by T.Takeda,    2011/06/01
|                                Y.Sato       2011/06/01
|                                K.Goto,      2010/01/12
|                                K.Matsubara, 2010/06/01
|
|   Contact address : IIS, The University of Tokyo CISS
|
 ----------------------------------------------------------
*/
namespace pmw
{
#ifndef EQNTERM_H_
#define EQNTERM_H_
class CEqnTerm
{
public:
	CEqnTerm() {}
	virtual ~CEqnTerm() {}
	void set(uiint meshID, uiint nodeID, uiint dof, double coef)
	{
		mMeshID = meshID;
		mNodeID = nodeID;
		mDOF = dof;
		mCoef = coef;
	}
	const uiint &meshID() const { return mMeshID; }
	const uiint &nodeID() const { return mNodeID; }
	const uiint &dof() const { return mDOF; }
	const double &coef() const { return mCoef; }
private:
	uiint mMeshID;
	uiint mNodeID;
	uiint mDOF;
	double mCoef;
};
#endif /* EQNTERM_H_ */
}
