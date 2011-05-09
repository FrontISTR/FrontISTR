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
	void set(int meshID, int nodeID, int dof, double coef)
	{
		mMeshID = meshID;
		mNodeID = nodeID;
		mDOF = dof;
		mCoef = coef;
	}
	const int &meshID() const { return mMeshID; }
	const int &nodeID() const { return mNodeID; }
	const int &dof() const { return mDOF; }
	const double &coef() const { return mCoef; }

private:
	int mMeshID;
	int mNodeID;
	int mDOF;
	double mCoef;
};
#endif /* EQNTERM_H_ */
}


