/*
 * EqnTerm.h
 *
 *  Created on: Nov 11, 2009
 *      Author: goto
 *
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


